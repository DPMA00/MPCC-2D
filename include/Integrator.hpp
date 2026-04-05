#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "Eigen/Dense"
#include "Model.hpp"

class Integrator
{
public:
    Integrator() = default;
    virtual ~Integrator() = default;

    virtual Eigen::VectorXd step(const Model& model, const Eigen::VectorXd& x, const Eigen::VectorXd& u, double Ts) const = 0;
    virtual void linearize_step(const Model& model, const Eigen::VectorXd& x, const Eigen::VectorXd&u,
                                double Ts, Eigen::MatrixXd& J_dyn,  Eigen::MatrixXd& J_func) const =0;
};

#endif