#include "MPCC.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <algorithm>

void write_paths_csv(const std::string& filename,
                     const std::vector<double>& xref,
                     const std::vector<double>& yref,
                     const std::vector<double>& xtraj,
                     const std::vector<double>& ytraj,
                     const std::vector<double>& velocity)
{
    std::ofstream file(filename);
    file << "xref,yref,xtraj,ytraj,vel\n";

    size_t nrows = std::max(xref.size(), xtraj.size());
    for (size_t i = 0; i < nrows; ++i) {
        if (i < xref.size()) file << xref[i];
        file << ",";

        if (i < yref.size()) file << yref[i];
        file << ",";

        if (i < xtraj.size()) file << xtraj[i];
        file << ",";

        if (i < ytraj.size()) file << ytraj[i];
        file << ",";

        if (i < velocity.size()) file << velocity[i];

        file << "\n";
    }
}


waypoints load_track_csv(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + filename);

    std::string line;

    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '#')
            break;
    }

    waypoints pts;


    auto parse_line = [&](const std::string& l) {
        std::istringstream ss(l);
        std::string token;
        double x, y;

        std::getline(ss, token, ','); x = std::stod(token);
        std::getline(ss, token, ','); y = std::stod(token);

        pts.x.push_back(x);
        pts.y.push_back(y);
    };

    if (!line.empty())
        parse_line(line);

    while (std::getline(file, line)) {
        if (!line.empty())
            parse_line(line);
    }

    return pts;
}

int main()
{
    ParametricSpline spline(T2_NATURAL_BOUNDARY_SPLINE);    
    MPCC mpcc(20, 0.1, spline);

    mpcc.configure_dynamics(SECOND_ORDER_MODEL, EXPL_EULER);
    mpcc.config_projection(NEWTON_STEP, 20, 1.0e-6);
    mpcc.config_solver_settings(1, 1e-3, 150, PRINT_LEVEL_SIMPLE);


    waypoints points = load_track_csv("/path/to/track.csv");
    
    mpcc.update_path(points);

    Eigen::VectorXd x0(7);
    x0 << points.x[0], points.y[0], M_PI/2, 0, 0, 0, 0;
    Eigen::VectorXd Q(2);
    Q<< 250, 50;

    Eigen::VectorXd R(4);
    R<<1,1,1,10;

    mpcc.set_weigths(Q, R);


    Eigen::VectorXd lbx(7);
    Eigen::VectorXd ubx(7);
    Eigen::VectorXd lbu(4);
    Eigen::VectorXd ubu(4);

    lbx << -1e6, -1e6, -M_PI,   0.0,  -4.0,  -0.8,   0.0;
    ubx <<  1e6,  1e6,  M_PI,  50.0,   4.0,   0.8,   1e6;

    lbu << -6.0, -5.0, -1.0,   0.0;
    ubu <<  3.0,  5.0,  1.0,  50.0;

    mpcc.set_constraints(lbx, ubx, lbu, ubu);
    int steps = 1000;//;
    std::vector<double> traj_x;
    std::vector<double> traj_y;
    std::vector<double> vel;

    traj_x.push_back(x0(0));
    traj_y.push_back(x0(1));

    for (int i = 0; i<steps; ++i)
    {
        mpcc.solve(x0);
        auto u = mpcc.get_controls();
        x0 = mpcc.simstep(x0, u);
        traj_x.push_back(x0(0));
        traj_y.push_back(x0(1));
        vel.push_back(x0(3));
    }
    double duplicate = vel.back();
    vel.push_back(duplicate);
    write_paths_csv("ref_and_traj.csv", points.x, points.y, traj_x, traj_y,vel);
    return 0;

}