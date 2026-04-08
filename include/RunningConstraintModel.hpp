#ifndef _RUNNING_CONSTRAINT_MODEL_H_
#define _RUNNING_CONSTRAINT_MODEL_H_

#include "Eigen/Dense"


using RowMajorMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class RunningConstraintModel
{

public:
    virtual ~RunningConstraintModel() = default;
    virtual int NrConstraints() const =0;

    virtual void add_constraint(RowMajorMat& A, Eigen::VectorXd& lbA, Eigen::VectorXd& ubA, int row, int colx, int colu, const Eigen::VectorXd& state, const Eigen::VectorXd& control) const = 0;
};

#endif
