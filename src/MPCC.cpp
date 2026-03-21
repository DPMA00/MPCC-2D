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
    J_R = Eigen::MatrixXd::Zero(2+NU,NVAR);
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

}


void MPCC::solve(const Eigen::VectorXd& x0)
{
    if (!current_path)
        throw std::runtime_error("build_nlprob: no path set");
    
    

    XY = x0.head<2>();
    s_prev = x0(NX-1);
    s = current_path->local_search(s_prev, XY);
    auto x_k = x0;
    x_k(NX-1) = s;

    X.block(0,0, NX,1) = x_k;

    auto x_k= x0;
    for (int k =0; k<N; ++k)
    {
        double x = x_k(0);
        double y = x_k(1);
        double theta = x_k(2);
        double s = x_k(3);
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

        double dcontour_s = cosphi*dphis - dxs*sinphi - xr*cosphi*dphis + sinphi*dphis + dys*cosphi - yr*sinphi*dphis;

        double dlag_s = cosphi*dphis - dxs*sinphi - xr*cosphi*dphis + sinphi*dphis + dys*cosphi - yr*sinphi*dphis;

    }
    

}