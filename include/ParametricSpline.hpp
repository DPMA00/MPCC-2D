#ifndef _PARAMETRIC_SPLINE_H_
#define _PARAMETRIC_SPLINE_H_

#include "Spline1D.hpp"


enum ProjMethod{
    NEWTON_STEP = 0,
    KDTREE = 1
};

struct NewtonConfig{
    double tolerance;
    int max_iter;
};

struct KDTreeConfig{
    double distance_upper_bound;
    double eps;
};


struct PathDat{
    double x;
    double y;
    double phi;
    double dxs;
    double dys;
    double dphis;
};

class ParametricSpline
{

private:
    std::vector<double> compute_arc_lengths(const waypoints& points);
    double newton_search(bool global,const double initial_guess, const Eigen::Vector2d& point);


    ProjMethod method;

    std::vector<double> s_;
    Spline1D spline_sx_;
    Spline1D spline_sy_;
    SplineType type_;

    NewtonConfig newton_config_;
    KDTreeConfig kdtree_config_;

public:
    ParametricSpline(SplineType splinetype);
    ~ParametricSpline();
    
    void set_proj_method(ProjMethod& method);
    void update_path(const waypoints& points);
    std::vector<double> get_arc_lengths() const;
    PathDat evalf_diff(double s);
    double local_search(const double initial_guess, const Eigen::Vector2d& point);
    void configure_newton(const NewtonConfig& config);
    void configure_kdtree(const KDTreeConfig& config);
};

#endif