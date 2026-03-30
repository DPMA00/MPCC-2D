#include "MPCC.hpp"
#include <iostream>


#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <algorithm>

void write_paths_csv(const std::string& filename,
                     const std::vector<double>& xref,
                     const std::vector<double>& yref,
                     const std::vector<double>& xtraj,
                     const std::vector<double>& ytraj)
{
    std::ofstream file(filename);

    file << "xref,yref,xtraj,ytraj\n";

    size_t nrows = std::max(xref.size(), xtraj.size());
    for (size_t i = 0; i < nrows; ++i) {
        if (i < xref.size()) file << xref[i];
        file << ",";

        if (i < yref.size()) file << yref[i];
        file << ",";

        if (i < xtraj.size()) file << xtraj[i];
        file << ",";

        if (i < ytraj.size()) file << ytraj[i];

        file << "\n";
    }
}

int main()
{
    ParametricSpline spline(T2_NATURAL_BOUNDARY_SPLINE);    
    MPCC mpcc(20, 0.1, spline);

    mpcc.configure_dynamics(DIFFDRIVE, EXPL_EULER);
    mpcc.config_projection(NEWTON_STEP, 20, 1.0e-6);
    mpcc.config_solver_settings(5, 1e-4, 100, 1e-6, PRINT_LEVEL_NONE);

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
    x0 << x[0], y[0], 0, 0;
    Eigen::VectorXd Q(2);
    Q<< 1000, 1000;

    Eigen::VectorXd R(3);
    R<<10,10,100;

    mpcc.set_weigths(Q, R);
    int steps = 100;//;
    std::vector<double> traj_x;
    std::vector<double> traj_y;

    traj_x.push_back(x0(0));
    traj_y.push_back(x0(1));

    for (int i = 0; i<steps; ++i)
    {
        mpcc.solve(x0);
        auto u = mpcc.get_solution();
        x0 = mpcc.simstep(x0, u);
        traj_x.push_back(x0(0));
        traj_y.push_back(x0(1));
    }

    write_paths_csv("ref_and_traj.csv", x, y, traj_x, traj_y);
    return 0;

}