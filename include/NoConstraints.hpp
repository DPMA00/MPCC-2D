#ifndef _NO_CONSTRAINTS_H_
#define _NO_CONSTRAINTS_H_

#include "RunningConstraintModel.hpp"


class NoConstraints : public RunningConstraintModel
{
private:
    int nrConstraints;
    
public:
    ~NoConstraints() override = default;
    int NrConstraints() const override {return 0;}
    void add_constraint(RowMajorMat& A, Eigen::VectorXd& lbA, Eigen::VectorXd& ubA, int row, int colx, int colu, const Eigen::VectorXd& state, const Eigen::VectorXd& control) const override;

};

#endif