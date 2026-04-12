#include "MPCC.hpp"

MPCC::MPCC(int N, double Ts, ParametricSpline& spline)
: N(N), Ts(Ts), current_path(&spline)
{
}

MPCC::~MPCC(){}



void MPCC::config_projection(ProjMethod proj, int max_iter, double tolerance)
{
    current_path->set_proj_method(proj);
    NewtonConfig config = {tolerance, max_iter};
    current_path->configure_newton(config);
}



void MPCC::config_projection(ProjMethod proj, double eps, double distance_upperbound)
{
    current_path->set_proj_method(proj);
    KDTreeConfig config = {distance_upperbound, eps};
    current_path->configure_kdtree(config);
}


void MPCC::config_solver_settings(int max_iter, double tol, int QP_max_iter, PrintLevel printsetting)
{
    sol_settings.max_iter = max_iter;
    sol_settings.tolerance = tol;
    sol_settings.QP_max_iter = QP_max_iter;
    sol_settings.printlevel = printsetting;

    

    qpOASES::Options options;
    options.setToMPC();
    options.enableEqualities = qpOASES::BT_TRUE;

    if (printsetting == PRINT_LEVEL_NONE)
    {
        options.printLevel = qpOASES::PL_NONE;
    }
    else if (printsetting==PRINT_LEVEL_SIMPLE)
    {
        options.printLevel = qpOASES::PL_LOW;
    }
    else if (printsetting==PRINT_LEVEL_DETAILED)
    {
        options.printLevel = qpOASES::PL_HIGH;
    }    

    qpprob.setOptions(options);
}

void MPCC::configure_dynamics(DynModel model, Integ integrator)
{
    switch (model)
    {
        case(DIFFDRIVE):
        {
            this->model = std::make_unique<DiffDriveModel>();
            runningcost = std::make_unique<DiffDriveCost>();
            constraints = std::make_unique<NoConstraints>();
            break;
        }
        case(EXTENDED_BICYCLE_MODEL):
        {
            this->model = std::make_unique<ExtendedBicycleModel>();
            runningcost = std::make_unique<ExtendedBicycleCost>();
            constraints = std::make_unique<NoConstraints>();
            break;//
        }
        case(SECOND_ORDER_MODEL):
        {   
            this->model = std::make_unique<SecondOrderModel>();
            runningcost = std::make_unique<SecondOrderModelCost>();
            constraints = std::make_unique<SecondOrderModelConstraints>();
            break;
        }
    }
    
    switch (integrator)
    {
        case(EXPL_EULER):
        {
            this->integrator = std::make_unique<ForwardEuler>();
            break;
        }
        case(EXPL_RK4):
        {
            //this->integrator = std::make_unique<RK4>();
            break;        
        }
    }
 
    NSLACK = 0;
    stage_const = constraints->NrConstraints(); // per stage
    NSLACK+=this->model->nS();
    W_S = this->model->W_S();


    NX = this->model->nx();
    NU = this->model->nu();
    NCON = NX*(N+1) + stage_const*N; // add equality constraints from multiple shooting and stage constraints for full prediction horizon
    NVAR = NX*(N+1) + NU*N; // multiple shooting without slack variables
    NVAR += NSLACK*N; // add slack variables;

    qpprob = qpOASES::SQProblem(NVAR,NCON);

    H =  RowMajorMat::Zero(NVAR,NVAR);
    g = Eigen::VectorXd::Zero(NVAR);
    A = RowMajorMat::Zero(NCON,NVAR);
    lbA = Eigen::VectorXd(NCON);
    ubA = Eigen::VectorXd(NCON);
    const double INF = 1e10;
    

    r_size = 2+NU; // contour error, lag error + controls

    W = Eigen::VectorXd(r_size);

    J_func = Eigen::MatrixXd::Zero(NX,NX+NU);
    J_dyn= Eigen::MatrixXd::Zero(NX,NX+NU);
    dyn_dev = Eigen::VectorXd::Zero(NX);

    r_z = Eigen::VectorXd::Zero(r_size);
    j_r = Eigen::MatrixXd::Zero(r_size,NX+NU);
    X = Eigen::MatrixXd::Zero(NX,N+1); // X linearization trajectory 
    U = Eigen::MatrixXd::Zero(NU,N); // U linearization trajectory
    dZ = Eigen::VectorXd::Zero(NVAR);
    NX_Identity = Eigen::MatrixXd::Identity(NX,NX);
    NU_Identity = Eigen::MatrixXd::Identity(NU,NU);


    S = Eigen::MatrixXd::Zero(NSLACK,N);


    lbz = Eigen::VectorXd::Constant(NVAR, -INF);
    ubz = Eigen::VectorXd::Constant(NVAR, INF);

    lbx = Eigen::VectorXd::Zero(NX);
    ubx = Eigen::VectorXd::Zero(NX);
    lbu = Eigen::VectorXd::Zero(NU);
    ubu = Eigen::VectorXd::Zero(NU);
    
    lbS = Eigen::VectorXd::Zero(NSLACK);
    ubS = Eigen::VectorXd::Constant(NSLACK,40); // will be modified later...
}


void MPCC::set_weigths(const Eigen::VectorXd& Qe, const Eigen::VectorXd& Ru)
{
    if (W.size() == Qe.size() + Ru.size())
    {
        W << Qe, Ru;
        // no terminal conditions added for now
    }
    else
    {
        throw std::runtime_error("set_weigths: Weight vector size mismatch");
    }
    
}


void MPCC::set_constraints(const Eigen::VectorXd& lbx, const Eigen::VectorXd&ubx,
                        const Eigen::VectorXd& lbu, const Eigen::VectorXd&ubu)
{
    this->lbx = lbx;
    this->ubx = ubx;
    this->lbu = lbu;
    this->ubu = ubu;

}



Eigen::VectorXd MPCC::RK4_step(const Eigen::VectorXd& x,const Eigen::VectorXd& u)
{
    Eigen::VectorXd k1 = model->dynamics(x, u);
    Eigen::VectorXd k2 = model->dynamics(x + 0.5*Ts * k1, u);
    Eigen::VectorXd k3 = model->dynamics(x + 0.5*Ts * k2, u);
    Eigen::VectorXd k4 = model->dynamics(x + Ts * k3, u);
    return x + Ts/6 *(k1 + 2*k2 +2*k3 + k4);
    

}


void MPCC::update_path(const waypoints& points)
{
    current_path->update_path(points);
}

void MPCC::warmstart(const Eigen::VectorXd& x0)
{
    // Shift the state and controls initial guess by removing the first entry nd duplicating the last
    
    Eigen::MatrixXd Xnew(NX,N+1);
    Xnew.leftCols(N) = X.rightCols(N);
    Xnew.col(N) = Xnew.col(N-1);
    X = Xnew;
    X.col(0) = x0;

    // Projection of X,Y onto track to retrieve closest path parameter s
    XY = x0.head<2>();
    s_prev = x0(NX-1);
    X(NX-1,0) = current_path->local_search(s_prev, XY); 
    
    


    Eigen::MatrixXd Unew(NU,N);
    Unew.leftCols(N-1) = U.rightCols(N-1);
    Unew.col(N-1) = Unew.col(N-2);
    U = Unew;
    

    Eigen::MatrixXd S_new(NSLACK,N);
    S_new.leftCols(N-1) = S.rightCols(N-1);
    S_new.col(N-1) = S_new.col(N-2);
    S = S_new;
    /*
    X.col(0) = x0;
    XY = x0.head<2>();
    s_prev = x0(NX-1);
    X(NX-1,0) = current_path->local_search(s_prev, XY); 
    
    for (int k = 0; k<N;++k)
    {
        X.col(k+1) = integrator->step(*model, X.col(k), U.col(k),Ts);
    }
    */
}




Eigen::VectorXd MPCC::get_controls()
{
    return U.col(0);
}

FullSolution MPCC::get_full_solution()
{
    FullSolution sol;
    sol.X = X;
    sol.U = U;
    return sol;
}


int MPCC::get_idx_x(int k)
{
    return k*(NX+NU);
}

int MPCC::get_idx_u(int k)
{
    return k*(NX+NU) + NX;
}

int MPCC::get_idx_S(int k)
{
    return NX*(N+1) + NU*N + NSLACK*k;
}

void MPCC::SQP_step(const Eigen::VectorXd& dz)
{
    double alpha= 1.0; //0.5, 0.8,1 works well
    for (int k=0; k<N; ++k)
    {
        X.col(k) += alpha*dz.segment(get_idx_x(k), NX);
        U.col(k) += alpha*dz.segment(get_idx_u(k), NU);
        S.col(k) += alpha*dz.segment(get_idx_S(k), NSLACK);
    }
    X.col(N) += alpha*dz.segment(get_idx_x(N), NX);
}

void MPCC::setupQP()
{
    const auto Wd = W.asDiagonal();
    const auto W_S_d = W_S.asDiagonal();
    int const_begin = NX*(N+1);
    for (int k =0; k<N; ++k)
        {
            // Linearization at z_k(bar) = [X_k, U_k]^T of the cost function + additional slack variables 
            auto state = X.col(k);
            auto control = U.col(k);
            auto slacks = S.col(k);

            // update jacobian and residuals and return pathdata (useful for curvature aware constraints)
            dat_s = runningcost->residual_and_jacobian(state, control, *current_path, r_z, j_r, dat_s);
            
            // Slacks Jacobian is zero in all entries except for the slacks at the current iterate k => identity matrix (NSLACK,NSLACK)

            // QP Cost
            g.segment(get_idx_x(k), NX+NU) += j_r.transpose() * Wd * r_z;
            H.block(get_idx_x(k), get_idx_x(k), NX+NU,NX+NU) += j_r.transpose() * Wd * j_r;
            g.segment(get_idx_S(k), NSLACK) += W_S_d * slacks;
            H.block(get_idx_S(k), get_idx_S(k), NSLACK, NSLACK) += W_S_d; 

            // multiple shooting constraints
            int mps_row = NX*(k+1);

            auto next_state = integrator->step(*model, state, control, Ts);
            dyn_dev = X.col(k+1) - next_state;
            integrator->linearize_step(*model, state, control, Ts, J_dyn, J_func);

            A.block(mps_row, (NX+NU)*k, NX, NX+NU) = J_func;
            A.block(mps_row, (NX+NU)*(k+1), NX,NX) = -NX_Identity;
            lbA.segment(mps_row, NX) = dyn_dev;
            ubA.segment(mps_row, NX) = dyn_dev;

            lbz.segment(get_idx_x(k),NX) = lbx - state;
            lbz.segment(get_idx_u(k),NU) = lbu - control;
            lbz.segment(get_idx_S(k),NSLACK) = lbS - slacks;

            ubz.segment(get_idx_x(k),NX) = ubx - state;
            ubz.segment(get_idx_u(k),NU) = ubu - control;
            ubz.segment(get_idx_S(k),NSLACK) = ubS - slacks;

            int row = const_begin + k*stage_const;

            constraints->add_constraint(A, lbA, ubA,row, get_idx_x(k),get_idx_u(k), get_idx_S(k),state,control,slacks, dat_s);
            j_r.setZero(); // reset the cost jacobians

            
        }
    lbz.segment(get_idx_x(N),NX) = lbx - X.col(N);
    ubz.segment(get_idx_x(N),NX) = ubx - X.col(N);
    H.diagonal().array() +=1e-6;
}

void MPCC::solve(const Eigen::VectorXd& x0)
{
    auto start = std::chrono::high_resolution_clock::now();
    if (!current_path)
        throw std::runtime_error("solve: no path set");
    
    warmstart(x0);
    bool init = false;
    qpOASES::returnValue ret;
    solved= false;
    //SQP Loop
    for (int i=0; i<sol_settings.max_iter; ++i)
    {
        
        H.setZero();
        g.setZero();
        A.setZero();
        lbA.setZero();
        ubA.setZero();

        // State at k=0 equality constraint
        A.block(0,0,NX,NX) = NX_Identity;
        lbA.segment(0,NX).setZero();
        ubA.segment(0,NX).setZero();
        int NWSR = sol_settings.QP_max_iter;
        
        setupQP();
        std::cout << "H diag min/max: " << H.diagonal().minCoeff() << " " << H.diagonal().maxCoeff() << "\n";
        if (!init)
        {
            ret = qpprob.init(H.data(),
                        g.data(),
                        A.data(),
                        lbz.data(),
                        ubz.data(),
                        lbA.data(),
                        ubA.data(), NWSR);
            init = (ret == qpOASES::SUCCESSFUL_RETURN);
        }

        else
        {
            ret = qpprob.hotstart(H.data(),
                        g.data(),
                        A.data(),
                        lbz.data(),
                        ubz.data(),
                        lbA.data(),
                        ubA.data(), NWSR);
            init = (ret == qpOASES::SUCCESSFUL_RETURN);
        }
        if (ret !=qpOASES::SUCCESSFUL_RETURN)
        {
            std::cout << "QP solver failed to find solution." << std::endl;
            break;
        }

        qpprob.getPrimalSolution(dZ.data());
        //std::cout << dZ.norm() << std::endl;
        if (dZ.norm()<sol_settings.tolerance)
        {
            std::cout << "Solution found!" << std::endl;
            solved = true;
            break;
        }
        SQP_step(dZ);

    }
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout << duration.count() << std::endl;
}


Eigen::VectorXd MPCC::simstep(const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
    return RK4_step(x, u);
}