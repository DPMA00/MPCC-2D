#ifndef _DYNAMIC_MODEL_H_
#define _DYNAMIC_MODEL_H_

#include "Eigen/Dense"


class Model
{
private:

public:
    virtual int nx() const = 0;
    virtual int nu() const = 0;
    virtual const Eigen::MatrixXd& I_NX() const = 0;
    virtual const Eigen::MatrixXd& I_NU() const = 0;


    virtual Eigen::VectorXd dynamics(const Eigen::VectorXd& x, const Eigen::VectorXd&u) const = 0;
    virtual void linearize_dynamics(const Eigen::VectorXd&x, const Eigen::VectorXd&u, Eigen::MatrixXd& J_dyn) const = 0;


    Model() = default;
    virtual ~Model() = default;
    
};

#endif