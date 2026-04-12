#ifndef _DIFFDRIVE_COST_H_
#define _DIFFDRIVE_COST_H_

#include "RunningCost.hpp"
#include "DiffDriveModel.hpp"

class DiffDriveCost : public RunningCost
{
public:

    ~DiffDriveCost() override = default;

    PathDat residual_and_jacobian(const Eigen::VectorXd& state,const Eigen::VectorXd& control, ParametricSpline& path,
        Eigen::VectorXd& res, Eigen::MatrixXd& j_r, PathDat& dat_s) const override;
};

#endif