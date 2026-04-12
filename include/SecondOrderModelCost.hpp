#ifndef _SECOND_ORDER_MODEL_COST_H_
#define _SECOND_ORDER_MODEL_COST_H_

#include "RunningCost.hpp"
#include "SecondOrderModel.hpp"

class SecondOrderModelCost : public RunningCost
{
public:

    ~SecondOrderModelCost() override = default;

    PathDat residual_and_jacobian(const Eigen::VectorXd& state,const Eigen::VectorXd& control, ParametricSpline& path,
        Eigen::VectorXd& res, Eigen::MatrixXd& j_r, PathDat& dat_s) const override;
};

#endif