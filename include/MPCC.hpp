#ifndef _MPCC_H_
#define _MPCC_H_

#include "ParametricSpline.hpp"
#include "qpOASES.hpp"
#include <optional>

enum DynModel{
    DIFFDRIVE = 0,
    BICYCLE_MODEL = 1
};

enum Integrator{
    EXPL_EULER =0,
    EXPL_RK4 = 1
};

class MPCC
{
private:
    ParametricSpline* current_path;
    DynModel model;
    Integrator integrator;
    PathDat dat_s; 

    

    int N;
    double Ts;
    int NX;
    int NU;
    int NCON;
    int NVAR;


    
    std::vector<double> H_dat;
    std::vector<double> g_dat;
    std::vector<double> A_dat;
    std::vector<double> lbA;
    std::vector<double> ubA;
    std::vector<double> lbX;
    std::vector<double> ubX;
    

    double s_prev;
    double s;
    Eigen::Vector2d XY;

    Eigen::MatrixXd X;
    Eigen::MatrixXd U;
    Eigen::MatrixXd J_R;

    Eigen::Vector4d diffdrive_dynamics(const Eigen::Vector4d& x,const Eigen::Vector3d&);
    Eigen::MatrixXd get_diffdrive_jacobian(const Eigen::VectorXd& x, const Eigen::VectorXd& u);

    Eigen::VectorXd RK4_step(const Eigen::VectorXd& X,const Eigen::VectorXd& U);


    void warmstart(const Eigen::VectorXd& x0);

    void getFullJacobian();

public:
    MPCC(int N, double Ts, ParametricSpline& spline);
    ~MPCC();
    void configure_dynamics(DynModel model = DIFFDRIVE, Integrator integrator = EXPL_RK4);
    void config_projection(ProjMethod proj, int max_iter, double tolerance);
    void config_projection(ProjMethod proj, double eps, double distance_upperbound);

    void solve(const Eigen::VectorXd& x0);
    void update_path(const waypoints& points);



};

#endif