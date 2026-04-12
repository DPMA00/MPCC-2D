#include "SecondOrderModel.hpp"


SecondOrderModel::SecondOrderModel()
: I_NX_(Eigen::MatrixXd::Identity(7,7)), I_NU_(Eigen::Matrix4d::Identity(4,4)), W_S_(Eigen::VectorXd(1))
{
    W_S_ << 1500;
}

Eigen::VectorXd SecondOrderModel::dynamics(const Eigen::VectorXd& x, const Eigen::VectorXd&u) const
{

    /*
        Dynamic model requires a constrained acceleration envelope based on experimental data.
        More details found in: https://arxiv.org/pdf/1703.01225 

    */
    double cospsi = std::cos(x(2));
    double sinpsi = std::sin(x(2));

    double vx = x(3);
    double vy = x(4);
    double vpsi = x(5);
    
    double xdot = vx*cospsi - vy*sinpsi;
    double ydot = vx*sinpsi + vy*cospsi; 
    double psidot = vpsi;

    double vxdot = u(0);
    double vydot = u(1);
    double vpsidot = u(2);
    double sdot = u(3);

    Eigen::VectorXd out(7);

    out << xdot, ydot, psidot, vxdot, vydot, vpsidot, sdot;

    return out;
}



void SecondOrderModel::linearize_dynamics(const Eigen::VectorXd&x, const Eigen::VectorXd&u, Eigen::MatrixXd& J_dyn) const
{
    J_dyn.setZero();
    
    double cospsi = std::cos(x(2));
    double sinpsi = std::sin(x(2));

    double vx = x(3);
    double vy = x(4);

    J_dyn(0,2) = -vx*sinpsi - vy*cospsi;
    J_dyn(0,3) = cospsi;
    J_dyn(0,4) = -sinpsi;

    J_dyn(1,2) = vx*cospsi - vy*sinpsi;
    J_dyn(1,3) = sinpsi;
    J_dyn(1,4) = cospsi;

    J_dyn(2,5) = 1.0;

    J_dyn(3,7) = 1.0;
    J_dyn(4,8) = 1.0;
    J_dyn(5,9) = 1.0;
    J_dyn(6,10) = 1.0;

}