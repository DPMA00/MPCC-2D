#include "DiffDriveModel.hpp"

DiffDriveModel::DiffDriveModel()
: I_NX_(Eigen::Matrix4d::Identity(4,4)), I_NU_(Eigen::Matrix3d::Identity(3,3)), W_S_(Eigen::VectorXd(0))
{}


Eigen::VectorXd DiffDriveModel::dynamics(const Eigen::VectorXd& x, const Eigen::VectorXd&u) const
{
    double xdot = u(0) * std::cos(x(2));
    double ydot = u(0) * std::sin(x(2));
    double thetadot = u(1);
    double sdot = u[2];
    Eigen::Vector4d out;
    out << xdot, ydot, thetadot,sdot;
    return out;
}

void DiffDriveModel::linearize_dynamics(const Eigen::VectorXd&x, const Eigen::VectorXd&u,
                           Eigen::MatrixXd& J_dyn) const
{
    J_dyn.setZero();

    double u1 = u(0);
    double x3 = x(2);
    double sintheta = std::sin(x3);
    double costheta = std::cos(x3);

    // A_c
    J_dyn(0,2) = -u1*sintheta;
    J_dyn(1,2) = u1*costheta;
    
    // B_c
    J_dyn(0,4) = costheta;
    J_dyn(1,4) = sintheta;
    J_dyn(2,5) = 1;
    J_dyn(3,6) = 1;
}