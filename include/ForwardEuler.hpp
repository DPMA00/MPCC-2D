#ifndef _FWD_EULER_INTEGRATOR_H_
#define _FWD_EULER_INTEGRATOR_H_

#include "Integrator.hpp"

class ForwardEuler : public Integrator
{
private:

public:
    ForwardEuler() = default;
    ~ForwardEuler() override = default;
    
    Eigen::VectorXd step(const Model& model, const Eigen::VectorXd& x, const Eigen::VectorXd& u, double Ts) const override;
    void linearize_step(const Model& model, const Eigen::VectorXd& x, const Eigen::VectorXd&u,
                                double Ts, Eigen::MatrixXd& J_dyn, Eigen::MatrixXd& J_func) const override;

};
#endif