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

    mpcc.configure_dynamics(DIFFDRIVE, EXPL_EULER);
    mpcc.config_projection(NEWTON_STEP, 20, 1.0e-6);
    mpcc.config_solver_settings(1, 1e-3, 50, PRINT_LEVEL_NONE);


    waypoints points = load_track_csv("/home/dpma/projects/mpcc_core/racetrack-database-master/tracks/Monza.csv");
    
    mpcc.update_path(points);

    Eigen::VectorXd x0(4);
    x0 << points.x[0], points.y[0], M_PI/2, 0;
    Eigen::VectorXd Q(2);
    Q<< 3300, 3000;

    Eigen::VectorXd R(3);
    R<<10,10,50;

    mpcc.set_weigths(Q, R);


    Eigen::VectorXd lbx(4);
    Eigen::VectorXd ubx(4);
    Eigen::VectorXd lbu(3);
    Eigen::VectorXd ubu(3);

    lbx << -INFINITY, -INFINITY, -INFINITY, -INFINITY;
    ubx << INFINITY, INFINITY, INFINITY, INFINITY;
    lbu << -20, -3, -70;
    ubu << 70,  3,  80;

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
        vel.push_back(u(0));
        if (x0(6) >= spline.get_path_length() -1.0)
        {        
            std::cout << "Lap time: " << (i+1)*0.1 << " seconds";
            break;
        }
    }
    double duplicate = vel.back();
    vel.push_back(duplicate);
    write_paths_csv("ref_and_traj.csv", points.x, points.y, traj_x, traj_y,vel);
    return 0;

}