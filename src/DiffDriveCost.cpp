#include "DiffDriveCost.hpp"

void DiffDriveCost::residual_and_jacobian(const Eigen::VectorXd&state,const Eigen::VectorXd& control, ParametricSpline& path, Eigen::VectorXd& res,
     Eigen::MatrixXd& j_r, PathDat& dat_s) const
{
    // Linearization at z_k(bar) = [X_k, U_k]^T of the cost function
        //z_k << state, control;
        double x = state(0);
        double y = state(1);
        double theta = state(2);
        double s = state(3);
        
        double u1 = control(0);
        double u2 = control(1);
        double uv = control(2); // virtual control

        dat_s = path.evalf_diff(s);
        double xr = dat_s.x;
        double yr = dat_s.y;
        double dxs = dat_s.dxs;
        double dys = dat_s.dys;
        double phi = dat_s.phi;
        double dphis = dat_s.dphis;

        double cosphi = std::cos(phi);
        double sinphi = std::sin(phi);

        // r_z(z_k) evaluation
        double e_c = (x-xr)*sinphi - (y-yr)*cosphi;
        double e_l = -(x-xr)*cosphi - (y-yr)*sinphi;
        double e_u1 = u1; // just minimize the controls for now (u-u_prev added later on)
        double e_u2 = u2; // ^==
        double e_uv = uv-80; // maximize progress (adjust target value according to performance behavior) 
        
        res << e_c, e_l, e_u1, e_u2, e_uv;


        // Jacobian j_r(z_k) entries 
        double dcontour_x = sinphi;
        double dcontour_y = -cosphi;
        double dcontour_s = x*cosphi*dphis - dxs*sinphi - xr*cosphi*dphis + y*sinphi*dphis + dys*cosphi - yr*sinphi*dphis;


        double dlag_x = -cosphi;
        double dlag_y = -sinphi;
        double dlag_s = x*sinphi*dphis + dxs*cosphi - xr*sinphi*dphis -y*cosphi*dphis + dys*sinphi + yr*cosphi*dphis;

        
        j_r(0,0) = dcontour_x;
        j_r(0,1) = dcontour_y;
        j_r(0,3) = dcontour_s;

        j_r(1,0) = dlag_x;
        j_r(1,1) = dlag_y;
        j_r(1,3) = dlag_s;
        j_r.block(2,4, 3,3) = Eigen::Matrix3d::Identity(3,3);
}