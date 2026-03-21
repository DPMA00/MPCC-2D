#ifndef _SPLINE1D_H_
#define _SPLINE1D_H_

#include "Eigen/Dense"
#include <vector>
#include <cmath>
#include <iostream>

struct waypoints{
    std::vector<double> x;
    std::vector<double> y;
};


struct segment{
    double x0;
    double a;
    double b;
    double c;
    double d;
};


enum SplineType{
    T2_NATURAL_BOUNDARY_SPLINE = 0,
    T1_BOUNDARY_SPLINE = 1,
    T2_BOUNDARY_SPLINE = 2
};

class Spline1D
{
private:
    double divided_diff(const int i_m,const int i,const int i_p);
    double get_mu(const int i_m,const int i,const int i_p);
    void solve_tridiag(const int size, const Eigen::VectorXd& mu, const Eigen::VectorXd& lambda,
                        const Eigen::VectorXd& d, Eigen::VectorXd& M);

    std::vector<segment> seg_;
    std::vector<double> X0;
    waypoints points_;
    size_t nrpoints;

public:

    Spline1D();
    void set_waypoints(waypoints points, SplineType type);
    void t2_natural_boundary_spline();

    int get_segment(double x);
    double evalf(int k, double x);
    double eval_derivative(int k, double x);
    double eval_second_derivative(int k, double x);

    ~Spline1D();
};



#endif
