#include "MPCC.hpp"

MPCC::MPCC(int N, double Ts, ParametricSpline& spline)
: N(N), Ts(Ts), current_path(&spline)
{

}

MPCC::~MPCC(){}



void MPCC::config_projection(ProjMethod proj, int max_iter, double tolerance)
{
    NewtonConfig config = {tolerance, max_iter};
    current_path->configure_newton(config);
}



void MPCC::config_projection(ProjMethod proj, double eps, double distance_upperbound)
{
    KDTreeConfig config = {distance_upperbound, eps};
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

    NCON = NX*N; // equality constraints from multiple shooting
    NVAR = NX*(N+1) + NU*N; // multiple shooting without slack variables
    H_dat.resize(NVAR*NVAR,0);
    g_dat.resize(NVAR,0);
    A_dat.resize(NCON*NVAR,0);
    lbA.resize(NCON,0);
    ubA.resize(NCON,0);
    lbX.resize(NVAR,0);
    ubX.resize(NVAR,0);

    J_dynamics = Eigen::MatrixXd(NX, NX+NU);
    J_R = Eigen::MatrixXd::Zero(2+NU,NVAR);
    j_r = Eigen::MatrixXd::Zero(2+NU,NX+NU);
    X = Eigen::MatrixXd::Zero(NX,N+1);
    U = Eigen::MatrixXd::Zero(NU,N);
}

Eigen::Vector4d MPCC::diffdrive_dynamics(const Eigen::Vector4d& X,const Eigen::Vector3d& U)
{
    double xdot = U(0) * std::cos(X(2));
    double ydot = U(0) * std::sin(X(2));
    double thetadot = U(1);
    double sdot = U[2];

    return Eigen::Vector4d {xdot, ydot, thetadot, sdot};
}

Eigen::VectorXd MPCC::RK4_step(const Eigen::VectorXd& x,const Eigen::VectorXd& u)
{
    if (model==DIFFDRIVE)
    {
        Eigen::VectorXd k1 = diffdrive_dynamics(x, u);
        Eigen::VectorXd k2 = diffdrive_dynamics(x + 0.5*Ts * k1, u);
        Eigen::VectorXd k3 = diffdrive_dynamics(x + 0.5*Ts * k2, u);
        Eigen::VectorXd k4 = diffdrive_dynamics(x + Ts * k3, u);
        return x + Ts/6 *(k1 + 2*k2 +2*k3 + k4);
    }
    throw std::runtime_error("RK4_step: unsupported dynamics model");
}


void MPCC::update_path(const waypoints& points)
{
    current_path->update_path(points);
}


void MPCC::warmstart(const Eigen::VectorXd& x0)
{
    //
}


void MPCC::linearize_model(const Eigen::VectorXd&X, const Eigen::VectorXd& U)
{
    if(model==DIFFDRIVE)
    {
        double u1 = U(0);
        double x3 = X(2);
        double sintheta = std::sin(x3);
        double costheta = std::cos(x3);
        Eigen::MatrixXd temp(2,3);
        temp << -u1*sintheta, 0, costheta,
                u1*costheta, 0, sintheta;

        Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
        J_dynamics.block(0,2,2,3) = temp;
        J_dynamics.block(2,5,2,2) = I;
    }
}


void MPCC::solve(const Eigen::VectorXd& x0)
{
    if (!current_path)
        throw std::runtime_error("build_nlprob: no path set");
    
    
    X.col(0) = x0;
    XY = x0.head<2>();
    s_prev = x0(NX-1);
    X(NX-1,0) = current_path->local_search(s_prev, XY);
    
    for (int k =0; k<N; ++k)
    {
        
        auto state = X.col(k);
        auto control = U.col(k);
        
        double x = state(0);
        double y = state(1);
        double theta = state(2);
        double s = state(3);
        

        dat_s = current_path->evalf_diff(s);
        double xr = dat_s.x;
        double yr = dat_s.y;
        double dxs = dat_s.dxs;
        double dys = dat_s.dys;
        double ddxs = dat_s.ddxs;
        double ddys = dat_s.ddys;
        double phi = dat_s.phi;
        double dphis = dat_s.dphis;

        double cosphi = std::cos(phi);
        double sinphi = std::sin(phi);

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

        j_r.block(2,NX, 3,3) = Eigen::MatrixXd::Identity(NU,NU);
        auto next_state = RK4_step(state, control);
        X.col(k+1) = next_state;

        J_R.block(0, (NX+NU)*k, 2+NU,NX+NU) = j_r;
        linearize_model(state, control);
        //
        
    }
    
    

}