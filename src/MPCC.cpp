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


void MPCC::config_solver_settings(int max_iter, double tol, int QP_max_iter, double QP_tol, PrintLevel printsetting)
{
    sol_settings.max_iter = max_iter;
    sol_settings.tolerance = tol;
    sol_settings.QP_max_iter = QP_max_iter;
    sol_settings.QP_tolerance = QP_tol;
    sol_settings.printlevel = printsetting;
}

void MPCC::configure_dynamics(DynModel model, Integrator integrator)
{
    this->model = model;
    this->integrator = integrator;

    if (model==DIFFDRIVE)
    {
        NX = 4;
        NU = 3;
    }
    else if (model==BICYCLE_MODEL)
    {
        NX = 6;
        NU = 3;
    }

    NCON = NX*(N+1); // equality constraints from multiple shooting
    NVAR = NX*(N+1) + NU*N; // multiple shooting without slack variables

    qpprob = qpOASES::SQProblem(NVAR,NCON);
    qpOASES::Options options;
    options.printLevel = qpOASES::PL_NONE;
    options.setToMPC();
    options.enableEqualities = qpOASES::BT_TRUE;
    //options.setToReliable();
    qpprob.setOptions(options);


    H =  RowMajorMat::Zero(NVAR,NVAR);
    g = Eigen::VectorXd::Zero(NVAR);
    A = RowMajorMat::Zero(NCON,NVAR);
    lbA = Eigen::VectorXd(NCON);
    ubA = Eigen::VectorXd(NCON);
    const double INF = 1e10;
    lbx = Eigen::VectorXd::Constant(NVAR, -INF);
    ubx = Eigen::VectorXd::Constant(NVAR, INF);

    r_size = 2+NU; // contour error, lag error + controls

    W = Eigen::VectorXd(r_size);

    J_func = Eigen::MatrixXd::Zero(NX,NX+NU);
    dyn_dev = Eigen::VectorXd::Zero(NX);

    r_z = Eigen::VectorXd::Zero(r_size);
    j_r = Eigen::MatrixXd::Zero(r_size,NX+NU);
    X = Eigen::MatrixXd::Zero(NX,N+1); // X linearization trajectory 
    U = Eigen::MatrixXd::Zero(NU,N); // U linearization trajectory
    dZ = Eigen::VectorXd::Zero(NVAR);
    NX_Identity = Eigen::MatrixXd::Identity(NX,NX);
    NU_Identity = Eigen::MatrixXd::Identity(NU,NU);
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


Eigen::Vector4d MPCC::diffdrive_dynamics(const Eigen::Vector4d& X,const Eigen::Vector3d& U)
{
    double xdot = U(0) * std::cos(X(2));
    double ydot = U(0) * std::sin(X(2));
    double thetadot = U(1);
    double sdot = U[2];
    Eigen::Vector4d out;
    out << xdot, ydot, thetadot,sdot;
    return out;
}

Eigen::VectorXd MPCC::Euler_step(const Eigen::VectorXd& x,const Eigen::VectorXd& u)
{
    if (model==DIFFDRIVE)
    {
        Eigen::Vector4d x_ = x.head<4>();
        Eigen::Vector3d u_ = u.head<3>();

        return x_ + Ts*diffdrive_dynamics(x_, u_);
    }
    throw std::runtime_error("Euler_step: unsupported dynamics model");
}


Eigen::VectorXd MPCC::RK4_step(const Eigen::VectorXd& x,const Eigen::VectorXd& u)
{
    if (model==DIFFDRIVE)
    {
        Eigen::Vector4d x_ = x.head<4>();
        Eigen::Vector3d u_ = u.head<3>();

        Eigen::Vector4d k1 = diffdrive_dynamics(x_, u_);
        Eigen::Vector4d k2 = diffdrive_dynamics(x_ + 0.5*Ts * k1, u_);
        Eigen::Vector4d k3 = diffdrive_dynamics(x_ + 0.5*Ts * k2, u_);
        Eigen::Vector4d k4 = diffdrive_dynamics(x_ + Ts * k3, u_);
        return x_ + Ts/6 *(k1 + 2*k2 +2*k3 + k4);
    }
    throw std::runtime_error("RK4_step: unsupported dynamics model");
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
    /*

    X.col(0) = x0;
    XY = x0.head<2>();
    s_prev = x0(NX-1);
    X(NX-1,0) = current_path->local_search(s_prev, XY); 
    
    for (int k = 0; k<N;++k)
    {
        X.col(k+1) = Euler_step(X.col(k), U.col(k));
    }*/
}


void MPCC::get_function_jacobian(const Eigen::VectorXd& X, const Eigen::VectorXd& U)
{
    if(integrator == EXPL_EULER)
    {
        if(model==DIFFDRIVE)
        {
            double u1 = U(0);
            double x3 = X(2);
            double sintheta = std::sin(x3);
            double costheta = std::cos(x3);

            J_func.setZero();

            // States: A_d = I + A_c * Ts
            J_func(0,0) = 1.0;
            J_func(1,1) = 1.0;
            J_func(2,2) = 1.0;
            J_func(3,3) = 1.0;

            J_func(0,2) = -u1*sintheta * Ts;
            J_func(1,2) = u1*costheta * Ts;

            // Controls: B_d = B_c*Ts;
            
            J_func(0,4) = costheta * Ts;
            J_func(1,4) = sintheta * Ts;
            J_func(2,5) = Ts;
            J_func(3,6) = Ts;

            return;
        }
    }
    throw std::runtime_error("get_function_jacobian: unsupported model/integrator");
}



Eigen::VectorXd MPCC::get_solution()
{
    return U.col(0);
}

int MPCC::get_idx_x(int k)
{
    return k*(NX+NU);
}

int MPCC::get_idx_u(int k)
{
    return k*(NX+NU) + NX;
}



void MPCC::SQP_step(const Eigen::VectorXd& dz)
{
    double alpha= 1; //0.5, 0.8,1 works well
    for (int k=0; k<N; ++k)
    {
        X.col(k) += alpha*dz.segment(get_idx_x(k), NX);
        U.col(k) += alpha*dz.segment(get_idx_u(k), NU);
    }
    X.col(N) += alpha*dz.segment(get_idx_x(N), NX);
}

void MPCC::setupQP()
{
    // Currently only supports diffdrive model with Euler 
    const auto Wd = W.asDiagonal();
    Eigen::VectorXd z_k(NX+NU);
    for (int k =0; k<N; ++k)
        {
            // Linearization at z_k(bar) = [X_k, U_k]^T of the cost function
            auto state = X.col(k);
            auto control = U.col(k);
            z_k << state, control;

            double x = state(0);
            double y = state(1);
            double theta = state(2);
            double s = state(3);
            
            double u1 = control(0);
            double u2 = control(1);
            double uv = control(2); // virtual control

            dat_s = current_path->evalf_diff(s);
            double xr = dat_s.x;
            double yr = dat_s.y;
            double dxs = dat_s.dxs;
            double dys = dat_s.dys;
            double phi = dat_s.phi;
            double dphis = dat_s.dphis;

            double cosphi = std::cos(phi);
            double sinphi = std::sin(phi);

            // r_z(z_k) evaluation
            double e_c = (x-xr)*sinphi - (y-yr)*cosphi;
            double e_l = -(x-xr)*cosphi - (y-yr)*sinphi;
            double e_u1 = u1; // just minimize the controls for now (u-u_prev added later on)
            double e_u2 = u2; // ^==
            double e_uv = uv-1; // maximize progress (adjust target value according to performance behavior) 
            
            r_z << e_c, e_l, e_u1, e_u2, e_uv;


            // Jacobian j_r(z_k) entries 
            double dcontour_x = sinphi;
            double dcontour_y = -cosphi;
            double dcontour_s = x*cosphi*dphis - dxs*sinphi - xr*cosphi*dphis + y*sinphi*dphis + dys*cosphi - yr*sinphi*dphis;


            double dlag_x = -cosphi;
            double dlag_y = -sinphi;
            double dlag_s = x*sinphi*dphis + dxs*cosphi - xr*sinphi*dphis -y*cosphi*dphis + dys*sinphi + yr*cosphi*dphis;

            
            j_r(0,0) = dcontour_x;
            j_r(0,1) = dcontour_y;
            j_r(0,3) = dcontour_s;

            j_r(1,0) = dlag_x;
            j_r(1,1) = dlag_y;
            j_r(1,3) = dlag_s;
            j_r.block(2,NX, NU,NU) = NU_Identity;

            g.segment(k*(NX+NU), NX+NU) += j_r.transpose() * Wd * r_z;

            H.block(k*(NX+NU), k*(NX+NU), NX+NU,NX+NU) += j_r.transpose() * Wd * j_r;
            

            auto next_state = Euler_step(state, control);


            // multiple shooting constraints
            dyn_dev = X.col(k+1) - next_state;
            get_function_jacobian(state,control);
            A.block(NX*(k+1), (NX+NU)*k, NX, NX+NU) = J_func;
            A.block(NX*(k+1), (NX+NU)*(k+1), NX,NX) = -NX_Identity;
            lbA.segment((k+1)*NX, NX) = dyn_dev;
            ubA.segment((k+1)*NX, NX) = lbA.segment((k+1)*NX,NX);

            j_r.setZero(); // reset the cost jacobians

            //
            
        }
    H.diagonal().array() +=1e-6;
}

void MPCC::solve(const Eigen::VectorXd& x0)
{
    auto start = std::chrono::high_resolution_clock::now();
    if (!current_path)
        throw std::runtime_error("solve: no path set");
    
    warmstart(x0);
    bool init = false;
    bool ret = false;
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
        if (!init)
        {
            ret = qpprob.init(H.data(),
                        g.data(),
                        A.data(),
                        lbx.data(),
                        ubx.data(),
                        lbA.data(),
                        ubA.data(), NWSR);
            init = (ret == qpOASES::SUCCESSFUL_RETURN);
        }

        else
        {
            ret = qpprob.hotstart(H.data(),
                        g.data(),
                        A.data(),
                        lbx.data(),
                        ubx.data(),
                        lbA.data(),
                        ubA.data(), NWSR);
            init = (ret == qpOASES::SUCCESSFUL_RETURN);
        }
        if (ret !=qpOASES::SUCCESSFUL_RETURN)
        {
            std::cout << "No solution found" << std::endl;
            break;
        }

        qpprob.getPrimalSolution(dZ.data());
        std::cout << dZ.norm() << std::endl;
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