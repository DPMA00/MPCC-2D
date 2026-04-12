#include "SecondOrderModelConstraints.hpp"

SecondOrderModelConstraints::SecondOrderModelConstraints() :  A_mat(Eigen::MatrixXd::Zero(6,4)), B_vec(Eigen::VectorXd::Zero(6))
{
    A_mat <<
    2.6,  1.0,  0.0, 0.0,
    2.6, -1.0,  0.0, 0.0,
    0.0,  1.1,  1.0, 0.0,
    0.0, -1.1, -1.0, 0.0,
    0.0, -0.57, 1.0, 0.0,
    0.0,  0.57,-1.0, 0.0;

    B_vec << 15.3, 15.3, 9.9, 9.9, 5.1, 5.1;
    
    alpha_ = 9.4;
    beta_ = 9.0;
}


void SecondOrderModelConstraints::add_constraint(RowMajorMat& A, Eigen::VectorXd& lbA, Eigen::VectorXd& ubA,
        int row, int colx, int colu, int colS, const Eigen::VectorXd& state, const Eigen::VectorXd& control,
        const Eigen::VectorXd& slack, PathDat& path) const
{
    double u1 = control(0);
    double u2 = control(1);
    double vx = state(3);
    double slack_vmax = slack(0);

    double kappa = std::abs(path.dphis);

    A.block(row, colu, 6,4) = A_mat;
    lbA.segment(row,6).setConstant(-1e10);
    ubA.segment(row,6) = B_vec-A_mat*control;

    double h_evalf_k = (u1*u1)/(alpha_*alpha_) + (u2*u2)/(beta_*beta_) -1.0;

    Eigen::RowVectorXd grad = Eigen::RowVectorXd::Zero(4);
    grad(0) =2.0*u1/(alpha_*alpha_);
    grad(1) = 2.0*u2/(beta_*beta_);

    A.block(row+6,colu, 1,4) = grad;
    lbA(row+6)= -1e10;
    ubA(row+6) = -h_evalf_k;


    double ax_min = -9.3 - 0.013*vx + 0.00072*vx*vx;
    double ax_max = 4.3 - 0.009*vx;

    A(row+7, colu+0) = 1.0;
    lbA(row+7) = -1e10;
    ubA(row+7) = ax_max-u1;

    A(row+8, colu+0) = -1.0;
    lbA(row+8) = -1e10;
    ubA(row+8) = -(ax_min-u1);


    double v_max_curv = std::sqrt(beta_ / std::max(kappa, 1e-4));
    A(row+9, colx+3) = 1.0;
    A(row+9, colS+0) = -1.0;
    lbA(row+9) = -1e10;
    ubA(row+9) = v_max_curv - vx + slack_vmax;
}