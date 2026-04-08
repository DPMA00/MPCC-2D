#ifndef _CONSTRAINED_2ND_ORDER_CONSTRAINTS_H_
#define _CONSTRAINED_2ND_ORDER_CONSTRAINTS_H_

#include "RunningConstraintModel.hpp"


class Constrained2ndOrderConstraints : public RunningConstraintModel
{
private:
    Eigen::MatrixXd A_mat;
    Eigen::VectorXd B_vec;
    double alpha_;
    double beta_;
    
public:
    ~Constrained2ndOrderConstraints() override = default;
    Constrained2ndOrderConstraints();
    int NrConstraints() const override {return 9;}
    void add_constraint(RowMajorMat& A, Eigen::VectorXd& lbA, Eigen::VectorXd& ubA, int row, int colx, int colu, const Eigen::VectorXd& state, const Eigen::VectorXd& control) const override;

};

#endif