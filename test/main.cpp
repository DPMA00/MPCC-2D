#include "MPCC.hpp"
#include <iostream>
#include <chrono>

int main()
{
    ParametricSpline spline(T2_NATURAL_BOUNDARY_SPLINE);    
    MPCC mpcc(50, 0.1, spline);

    mpcc.configure_dynamics(DIFFDRIVE, EXPL_RK4);
    mpcc.config_projection(NEWTON_STEP, 100, 1.0e-6);

    waypoints points;
    std::vector<double> x{0.0, 0.5, 0.7, 0.9, 1.5};
    std::vector<double> y{0.0, 0.3, 0.8, 1.3, 1.45};

    points.x = x;
    points.y = y;

    mpcc.update_path(points);

    Eigen::VectorXd x0(4);
    x0 << 0.46, 0.24, 0.0, 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    mpcc.solve(x0);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout << duration.count() << std::endl;
    
    return 0;

}