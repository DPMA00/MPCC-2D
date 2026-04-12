#include "NoConstraints.hpp"

void NoConstraints::add_constraint(RowMajorMat& A, Eigen::VectorXd& lbA, Eigen::VectorXd& ubA,
        int row, int colx, int colu, int colS, const Eigen::VectorXd& state, const Eigen::VectorXd& control,
        const Eigen::VectorXd& slack, PathDat& path) const
{
    
}