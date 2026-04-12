#ifndef _CONSTRAINED_2ND_ORDER_CONSTRAINTS_H_
#define _CONSTRAINED_2ND_ORDER_CONSTRAINTS_H_

#include "RunningConstraintModel.hpp"


class SecondOrderModelConstraints : public RunningConstraintModel
{
private:
    Eigen::MatrixXd A_mat;
    Eigen::VectorXd B_vec;
    double alpha_;
    double beta_;
    
public:
    ~SecondOrderModelConstraints() override = default;
    SecondOrderModelConstraints();
    int NrConstraints() const override {return 10;}
    void add_constraint(RowMajorMat& A, Eigen::VectorXd& lbA, Eigen::VectorXd& ubA,
        int row, int colx, int colu, int colS, const Eigen::VectorXd& state, const Eigen::VectorXd& control,
        const Eigen::VectorXd& slack, PathDat& path) const override;

};

#endif