#include "ParametricSpline.hpp"

ParametricSpline::ParametricSpline(SplineType splinetype)
: type_(splinetype)
{
    newton_config_.max_iter = 100;
    newton_config_.tolerance = 1e-6;
    kdtree_config_.distance_upper_bound = 0.005;
    kdtree_config_.eps = 0;
}

ParametricSpline::~ParametricSpline()
{}

void ParametricSpline::configure_newton(const NewtonConfig& config)
{
    newton_config_ = config;
}


std::vector<double> ParametricSpline::compute_arc_lengths(const waypoints& points)
{
    /*
    \brief 
    */
    const auto& x = points.x;
    const auto& y = points.y;

    size_t nrwaypoints = x.size();
    std::vector<double> s(nrwaypoints);

    s[0] = 0.0;

    for (size_t i=1; i<nrwaypoints; ++i)
    {
        double dx = x[i] - x[i-1];
        double dy = y[i] - y[i-1];
        s[i] = s[i-1] + std::sqrt(dx*dx + dy*dy);
    }
    return s;
}


void ParametricSpline::update_path(const waypoints& points)
{
    s_ = compute_arc_lengths(points);
    waypoints points_sx;
    waypoints points_sy;
    
    points_sx.x = s_;
    points_sx.y = points.x;
    
    points_sy.x = s_;
    points_sy.y = points.y;
    
   
    spline_sx_.set_waypoints(points_sx, type_);
    spline_sy_.set_waypoints(points_sy, type_);

}


std::vector<double> ParametricSpline::get_arc_lengths() const
{
    return s_;
}


PathDat ParametricSpline::evalf_diff(double s)
{
    int k = spline_sx_.get_segment(s);

    PathDat eval;
    eval.x = spline_sx_.evalf(k,s);
    eval.y = spline_sy_.evalf(k,s);
    eval.dxs = spline_sx_.eval_derivative(k,s);
    eval.dys = spline_sy_.eval_derivative(k,s);
    eval.ddxs = spline_sx_.eval_second_derivative(k,s);
    eval.ddys = spline_sy_.eval_second_derivative(k,s);
    eval.phi = std::atan2(eval.dys,eval.dxs);
    eval.dphis = (eval.dxs*eval.ddys - eval.dys*eval.ddxs) / (eval.dxs * eval.dxs + eval.dys*eval.dys +1e-9);
    return eval;
}


double ParametricSpline::local_search(const double initial_guess, const Eigen::Vector2d& point)
{
    switch(method)
    {
        case NEWTON_STEP:
        return newton_search(false, initial_guess, point);
        break;

        case KDTREE:
        return 2.0;
        break;
    }

    return 2.0;
}

double ParametricSpline::newton_search(bool global,const double initial_guess, const Eigen::Vector2d& point)
{
    /*
    \brief Performs newton iterations over 

    */
    
    const double x = point(0);
    const double y = point(1);
    
    double s = initial_guess;

    for (int i=0; i<newton_config_.max_iter; ++i)
    {
        int k = spline_sx_.get_segment(s);
        double Cx = spline_sx_.evalf(k, s);
        double Cy = spline_sy_.evalf(k, s);

        double dCx = spline_sx_.eval_derivative(k, s);
        double dCy = spline_sy_.eval_derivative(k, s);

        double ddCx = spline_sx_.eval_second_derivative(k, s);
        double ddCy = spline_sy_.eval_second_derivative(k, s);

        double dR = 2*((Cx -x)*dCx + (Cy -y)*dCy);
        double ddR =  2*(dCx*dCx + (Cx-x)*ddCx + dCy*dCy + (Cy-y)*ddCy);
        
        if (std::abs(ddR) < 1e-12) break;
        
        double step = dR/ddR;
        double s_new = s-step;

        s = s_new;

        if (std::abs(step) < newton_config_.tolerance) break;
    }

    return s;
}
