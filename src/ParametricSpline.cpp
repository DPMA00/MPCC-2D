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
    /*
        Defines configuration settings for the newton type projection method.
    */
    newton_config_ = config;
}


void ParametricSpline::configure_kdtree(const KDTreeConfig& config)
{
    /*
        Defines configuration settings for the newton type projection method.
    */
    kdtree_config_ = config;
}


void ParametricSpline::set_proj_method(ProjMethod& method)
{
    this->method = method;
}

std::vector<double> ParametricSpline::compute_arc_lengths(const waypoints& points)
{
    /*
        Computes the approximate arc lengths up until the last waypoint by calculating the euclidean norm
        between consecutive waypoints. 
        
        Each cumulative arclength is stored in the new parameter s.
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
    /*
        Updates parametric spline representation of the path from the given waypoints.
        
        The waypoint sequence is first converted into arc-length parameter s.
        The x and y coordinates are then re-parametrized with respect to s, and two cubic splines are constructed:
            -x(s)
            -y(s)
    */
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
    // Returns the parameter vector s.
    return s_;
}


PathDat ParametricSpline::evalf_diff(double s)
{
    /*
        Evaluates necessary function values, derivatives and second derivatives of the reference path at s.

        The spline segment where the value of s corresponds with is first determined to evaluate all PathDat entries:

            - X_ref(s)
            - Y_ref(s)
            - ∂X_ref(s)/∂s
            - ∂Y_ref(s)/∂s  
            - φ(s)
            - ∂²X_ref(s)/∂s²
            - ∂²Y_ref(s)/∂s²
            - ∂φ(s)/∂s
        
            These are then returned into Pathdat eval.        

    */
    int k = spline_sx_.get_segment(s);

    double ddxs = spline_sx_.eval_second_derivative(k,s);
    double ddys = spline_sy_.eval_second_derivative(k,s);

    PathDat eval;
    eval.x = spline_sx_.evalf(k,s);
    eval.y = spline_sy_.evalf(k,s);
    eval.dxs = spline_sx_.eval_derivative(k,s);
    eval.dys = spline_sy_.eval_derivative(k,s);
    eval.phi = std::atan2(eval.dys,eval.dxs);
    eval.dphis = (eval.dxs*ddys - eval.dys*ddxs) / (eval.dxs * eval.dxs + eval.dys*eval.dys +1e-9);
    return eval;
}


double ParametricSpline::local_search(const double initial_guess, const Eigen::Vector2d& point)
{
    /*
        The selected search method evaluates the nearest path parameter s to the given point
        using a local search around the initial guess for s.
    */
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
        Uses Newton's method to compute the path parameter s corresponding to a local projection of the query point
        on to the path.
        
        The method a small nonlinear optimization problem that seeks to minimize the squared
        Euclidean distance between the query point (x,y) and the path point C(s) = [C_x(s), C_y(s)]:

                                R(s) = (x-C_x(s))^2 + (y-C_y(s))^2

        A Newton iteration is applied to the first order optimality condition
                                        R'(s) = 0
        using the update:
        s_{i+1} = s_i - R'(s_i)/R''(s_i)

        The initial guess determines the starting point of the search and therefore can strongly influence the
        solution obtained. At each iteration the spline segment corresponding to the current value of s is used
        to evaluate the path and its derivatives.
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
        s = std::clamp(s_new, s_.front(), s_.back()); // dont accidentally jump off the segment
        if (std::abs(step) < newton_config_.tolerance) break;
    }

    return s;
}
