#ifndef _DIFFERENTIAL_DRIVE_MODEL_H_
#define _DIFFERENTIAL_DRIVE_MODEL_H_

#include "Model.hpp"

class DiffDriveModel : public Model
{
private:
    Eigen::MatrixXd I_NX_;
    Eigen::MatrixXd I_NU_;

public:
    int nx() const override {return 4;}
    int nu() const override {return 3;}
    const Eigen::MatrixXd& I_NX() const override {return I_NX_;}
    const Eigen::MatrixXd& I_NU() const override {return I_NU_;}

 
    Eigen::VectorXd dynamics(const Eigen::VectorXd& x, const Eigen::VectorXd&u) const override;
    void linearize_dynamics(const Eigen::VectorXd&x, const Eigen::VectorXd&u, Eigen::MatrixXd& J_dyn) const override;
    
    DiffDriveModel();
    ~DiffDriveModel() override = default;
};


#endif