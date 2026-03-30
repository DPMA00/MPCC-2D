#ifndef _MPCC_H_
#define _MPCC_H_

#include "ParametricSpline.hpp"
#include "qpOASES.hpp"
#include <optional>
#include <chrono>


using RowMajorMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

enum DynModel{
    DIFFDRIVE = 0,
    BICYCLE_MODEL = 1
};

enum Integrator{
    EXPL_EULER =0,
    EXPL_RK4 = 1
};

enum PrintLevel{
    PRINT_LEVEL_NONE = 0,
    PRINT_LEVEL_SIMPLE =1,
    PRINT_LEVEL_DETAILED =2
};

struct SolverSettings{
    int max_iter;
    double tolerance;
    int QP_max_iter;
    double QP_tolerance;
    PrintLevel printlevel;
};

class MPCC
{
private:
    ParametricSpline* current_path;
    DynModel model;
    Integrator integrator;
    PathDat dat_s; 
    SolverSettings sol_settings;
    qpOASES::QProblem qpprob;
    

    int N;
    double Ts;
    int NX;
    int NU;
    int NCON;
    int NVAR;
    int r_size;

    
    double s_prev;
    double s;
    Eigen::Vector2d XY;

    Eigen::MatrixXd X;
    Eigen::MatrixXd U;

    Eigen::VectorXd dZ;

    Eigen::MatrixXd NX_Identity;
    Eigen::MatrixXd NU_Identity;

    RowMajorMat H;
    Eigen::VectorXd g;
    RowMajorMat A;
    Eigen::VectorXd lbA;
    Eigen::VectorXd ubA;
    Eigen::VectorXd lbx;
    Eigen::VectorXd ubx;

    Eigen::VectorXd W;
    Eigen::VectorXd r_z;
    Eigen::MatrixXd j_r;

    Eigen::MatrixXd J_func;
    Eigen::VectorXd dyn_dev;

    Eigen::Vector4d diffdrive_dynamics(const Eigen::Vector4d& x,const Eigen::Vector3d&);
    Eigen::MatrixXd get_diffdrive_jacobian(const Eigen::VectorXd& x, const Eigen::VectorXd& u);

    Eigen::VectorXd RK4_step(const Eigen::VectorXd& X,const Eigen::VectorXd& U);
    Eigen::VectorXd Euler_step(const Eigen::VectorXd& X, const Eigen::VectorXd& U);
    int get_idx_x(int k);
    int get_idx_u(int k);

    void SQP_step(const Eigen::VectorXd& dz);
    void setupQP();
    void warmstart(const Eigen::VectorXd& x0);
    void get_function_jacobian(const Eigen::VectorXd&X, const Eigen::VectorXd& U);

public:
    MPCC(int N, double Ts, ParametricSpline& spline);
    ~MPCC();
    void configure_dynamics(DynModel model = DIFFDRIVE, Integrator integrator = EXPL_RK4);
    void config_projection(ProjMethod proj, int max_iter, double tolerance);
    void config_projection(ProjMethod proj, double eps, double distance_upperbound);
    void config_solver_settings(int max_iter, double tol, int QP_max_iter, double QP_tol, PrintLevel printsetting);
    void set_weigths(const Eigen::VectorXd& Qe, const Eigen::VectorXd& Ru);

    void solve(const Eigen::VectorXd& x0);
    void update_path(const waypoints& points);

    Eigen::VectorXd get_solution();
    Eigen::VectorXd simstep(const Eigen::VectorXd& x, const Eigen::VectorXd& u);


};

#endif