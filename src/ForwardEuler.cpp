#include "ForwardEuler.hpp"

Eigen::VectorXd ForwardEuler::step(const Model& model, const Eigen::VectorXd& x, const Eigen::VectorXd& u, double Ts) const
{
    return x + model.dynamics(x,u)*Ts;
}


void ForwardEuler::linearize_step(const Model& model, const Eigen::VectorXd& x, const Eigen::VectorXd&u,
                                double Ts,  Eigen::MatrixXd& J_dyn, Eigen::MatrixXd& J_func) const
{
    /*
        Get discrete time Jacobian from the continuous time Jacobian after Euler integration:

        xdot(t) = f(x,u), let z = [x;u] and J_dyn = df/dz = [A_c, B_c]
        
        Forward Euler discretization:
       
        f(z_k) = x_k + f_dyn(z_k) * Ts

        J(z_k) = [I 0] + J_dyn * Ts
    */
   
    model.linearize_dynamics(x,u,J_dyn);
    J_func.noalias() = Ts* J_dyn;
    J_func.leftCols(model.nx()).diagonal().array() += 1.0;
}