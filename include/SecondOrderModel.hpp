#ifndef _CONSTRAINED_2ND_ORDER_MODEL_H_
#define _CONSTRAINED_2ND_ORDER_MODEL_H_

#include "Model.hpp"

class SecondOrderModel : public Model
{
private:
    Eigen::MatrixXd I_NX_;
    Eigen::MatrixXd I_NU_;
    Eigen::VectorXd W_S_;
public:

    int nx() const override {return 7;}
    int nu() const override {return 4;}
    int nS() const override {return 1;}
    Eigen::VectorXd W_S() const override {return W_S_;}
    const Eigen::MatrixXd& I_NX() const override {return I_NX_;}
    const Eigen::MatrixXd& I_NU() const override {return I_NU_;}

 
    Eigen::VectorXd dynamics(const Eigen::VectorXd& x, const Eigen::VectorXd&u) const override;
    void linearize_dynamics(const Eigen::VectorXd&x, const Eigen::VectorXd&u, Eigen::MatrixXd& J_dyn) const override;
    
    SecondOrderModel();
    ~SecondOrderModel() override = default;
};

#endif
