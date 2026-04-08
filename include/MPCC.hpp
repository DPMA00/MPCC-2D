#ifndef _MPCC_H_
#define _MPCC_H_

#include "ParametricSpline.hpp"
#include "qpOASES.hpp"
#include "ForwardEuler.hpp"
#include "DiffDriveModel.hpp"
#include "DiffDriveCost.hpp"
#include "Constrained2ndOrderModel.hpp"
#include "SecondOrderModelCost.hpp"
#include "ExtendedBicycleModel.hpp"
#include "ExtendedBicycleCost.hpp"
#include "NoConstraints.hpp"
#include "Constrained2ndOrderConstraints.hpp"
#include <memory>
#include <chrono>


using RowMajorMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

enum DynModel{
    DIFFDRIVE = 0,
    EXTENDED_BICYCLE_MODEL = 1,
    SECOND_ORDER_MODEL = 2
};

enum Integ{
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
    PrintLevel printlevel;
};


struct FullSolution{
    Eigen::VectorXd X;
    Eigen::VectorXd U;
};

class MPCC
{
private:
    std::unique_ptr<Model> model;
    std::unique_ptr<Integrator> integrator;
    std::unique_ptr<RunningCost> runningcost;
    std::unique_ptr<RunningConstraintModel> constraints;
    ParametricSpline* current_path;


    PathDat dat_s; 
    SolverSettings sol_settings;
    qpOASES::SQProblem qpprob;
    

    int N;
    double Ts;
    int stage_const;
    int NX;
    int NU;
    int NCON;
    int NVAR;
    int r_size;
    bool solved;
    
    double s_prev;
    double s;
    Eigen::Vector2d XY;

    Eigen::MatrixXd X;
    Eigen::MatrixXd U;

    Eigen::VectorXd dZ;

    Eigen::MatrixXd NX_Identity;
    Eigen::MatrixXd NU_Identity;


    // QP Params
    RowMajorMat H;
    Eigen::VectorXd g;
    RowMajorMat A;
    Eigen::VectorXd lbA;
    Eigen::VectorXd ubA;
    Eigen::VectorXd lbz;
    Eigen::VectorXd ubz;

    Eigen::VectorXd lbx;
    Eigen::VectorXd ubx;
    Eigen::VectorXd lbu;
    Eigen::VectorXd ubu;

    Eigen::VectorXd W;
    Eigen::VectorXd r_z;
    Eigen::MatrixXd j_r;

    Eigen::MatrixXd J_func;
    Eigen::MatrixXd J_dyn;
    Eigen::VectorXd dyn_dev;


    Eigen::VectorXd RK4_step(const Eigen::VectorXd& X,const Eigen::VectorXd& U);

    int get_idx_x(int k);
    int get_idx_u(int k);

    void SQP_step(const Eigen::VectorXd& dz);
    void setupQP();
    void warmstart(const Eigen::VectorXd& x0);


public:
    MPCC(int N, double Ts, ParametricSpline& spline);
    ~MPCC();
    void configure_dynamics(DynModel model = DIFFDRIVE, Integ integrator = EXPL_RK4);
    void config_projection(ProjMethod proj, int max_iter, double tolerance);
    void config_projection(ProjMethod proj, double eps, double distance_upperbound);
    void config_solver_settings(int max_iter, double tol, int QP_max_iter, PrintLevel printsetting);
    
    void set_weigths(const Eigen::VectorXd& Qe, const Eigen::VectorXd& Ru);
    void set_constraints(const Eigen::VectorXd& lbx, const Eigen::VectorXd&ubx,
                        const Eigen::VectorXd& lbu, const Eigen::VectorXd&ubu);

    void solve(const Eigen::VectorXd& x0);
    void update_path(const waypoints& points);

    Eigen::VectorXd get_controls();
    FullSolution get_full_solution();
    Eigen::VectorXd simstep(const Eigen::VectorXd& x, const Eigen::VectorXd& u);


};

#endif