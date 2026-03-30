#include "MPCC.hpp"
#include <iostream>

int main()
{
    ParametricSpline spline(T2_NATURAL_BOUNDARY_SPLINE);    
    MPCC mpcc(20, 0.1, spline);

    mpcc.configure_dynamics(DIFFDRIVE, EXPL_EULER);
    mpcc.config_projection(NEWTON_STEP, 20, 1.0e-6);
    mpcc.config_solver_settings(1, 1e-6, 50, 1e-6, PRINT_LEVEL_NONE);

    waypoints points;
    std::vector<double> x,y;
    int nrpoints = 200;

    for (int i=0; i<nrpoints; ++i)
    {
        double t = 2.0 * M_PI * i /nrpoints;
        x.push_back(20.0*std::cos(t));
        y.push_back(10.0*std::sin(2*t));
    }

    points.x = x;
    points.y = y;

    mpcc.update_path(points);

    Eigen::VectorXd x0(4);
    x0 << 0.46, 0.24, 0.0, 0.0;
    Eigen::VectorXd Q(2);
    Q<< 100, 100;

    Eigen::VectorXd R(3);
    R<<10,10,10;

    mpcc.set_weigths(Q, R);
    int steps = 5;
    //auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i<steps; ++i)
    {
        mpcc.solve(x0);
        auto u = mpcc.get_solution();
        x0 = mpcc.simstep(x0, u);
    }
    return 0;

}