#ifndef _RUNNING_COST_H_
#define _RUNNING_COST_H_

#include "Eigen/Dense"
#include "ParametricSpline.hpp"


class RunningCost
{
public:
    virtual ~RunningCost() = default;

    virtual PathDat residual_and_jacobian(const Eigen::VectorXd& state,const Eigen::VectorXd& control,
                ParametricSpline& path, Eigen::VectorXd& res, Eigen::MatrixXd& j_r, PathDat& dat_s) const = 0;

};


#endif