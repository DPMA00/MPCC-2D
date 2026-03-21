#include "Spline1D.hpp"


Spline1D::Spline1D()
{}

Spline1D::~Spline1D(){}

void Spline1D::set_waypoints(waypoints points, SplineType type)
{
    /* 
        Moves points into class and initialize nr segments and nr of waypoints,
        then calls selected spline constructor
    */
    points_ = std::move(points);
    nrpoints = points_.x.size();
    seg_.resize(nrpoints-1);


    switch(type)
    {
        case T2_NATURAL_BOUNDARY_SPLINE:
            t2_natural_boundary_spline();
            break;

        case T1_BOUNDARY_SPLINE:
            break;
            
        case T2_BOUNDARY_SPLINE: 
            break;
    }
}





double Spline1D::divided_diff(const int i_m,const int i,const int i_p)
{
    // Calculates a three point divided difference
    const auto& x = points_.x;
    const auto& y = points_.y;

    double left_term = (y[i_p] - y[i]) / (x[i_p]-x[i]);
    double right_term = (y[i] - y[i_m]) / (x[i]-x[i_m]);
    return (left_term-right_term) / (x[i_p]-x[i_m]);
}


double Spline1D::get_mu(const int i_m,const int i,const int i_p)
{
    // Calculates mu_i 
    const auto& x = points_.x;
    double num = x[i] - x[i_m];
    double denom = - x[i_m] + x[i_p]; //  x[i] - x[i-1] + x[i+1] - x[i]

    return num/denom;
}


void Spline1D::solve_tridiag(const int size, const Eigen::VectorXd& mu, const Eigen::VectorXd& lambda,
                        const Eigen::VectorXd& d, Eigen::VectorXd& M)
{
    // Solves the tridiagonal linear system using Thomas' algorithm
    Eigen::VectorXd scratch(size);
    M.resize(size);

    scratch(0) = lambda(0) / 2.0;
    M(0) = d(0) / 2.0;
    
    for (int i=1; i<size ; ++i)
    {
        const double denom = 2.0 - mu(i)* scratch(i-1);
        scratch(i) = lambda(i) / denom;
        
        M(i) = (d(i) - mu(i)*M(i-1)) / denom;
    }

    for (int i = size-2; i>=0; --i)
    {
        M(i) -= scratch(i) * M(i+1);
    }
}

void Spline1D::t2_natural_boundary_spline()
{
    /* 
    Calculates the 3rd order polynomial coefficients (a_i, b_i, c_i,d_i) for a type 2 normal boundary spline
    for each segment:
            C_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3     for i = 0...n-2
            
            https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation
    */
    const auto& x = points_.x;
    const auto& y = points_.y;

    Eigen::VectorXd M(nrpoints);
    Eigen::VectorXd d = Eigen::VectorXd::Zero(nrpoints);

    Eigen::VectorXd mu = Eigen::VectorXd::Zero(nrpoints);
    Eigen::VectorXd lambda = Eigen::VectorXd::Zero(nrpoints);


    for (int i=1; i<nrpoints-1; ++i)
    {
        d(i) = 6*divided_diff(i-1, i, i+1);
        mu(i) = get_mu(i-1, i, i+1);
        lambda(i) = 1.0-mu(i);
    }

    solve_tridiag(nrpoints, mu, lambda, d, M);

    for (int i=0 ; i< nrpoints-1; ++i)
    {
        double h = x[i+1] - x[i];
        double a_val = y[i];
        double b_val = (y[i+1]-y[i])/h - h/6 * (2*M(i)+M(i+1));
        double c_val = M(i)/2;
        double d_val = (M(i+1) - M(i))/(6*h);

        seg_[i] = {x[i], a_val, b_val, c_val, d_val};
    }

}

int Spline1D::get_segment(double x)
{
    // Get the index of the segment wherein x lies
    auto it = std::lower_bound(points_.x.begin(), points_.x.end(),x);
    int k = std::distance(points_.x.begin(), it) -1;
    k = std::max(0, std::min(k, (int)seg_.size() -1));

    return k;
}

double Spline1D::evalf(int k, double x)
{
    // Evaluates the spline segment C_k at x

    const double& a = seg_[k].a;
    const double& b = seg_[k].b;
    const double& c = seg_[k].c;
    const double& d = seg_[k].d;

    double h = x-seg_[k].x0;

    return a + b*h + c*std::pow(h, 2) + d*std::pow(h,3);
    
}


double Spline1D::eval_derivative(int k, double x)
{
    // Evaluates the spline segment derivative C_k' at x
    const double& x0 = seg_[k].x0;
    const double& a = seg_[k].a;
    const double& b = seg_[k].b;
    const double& c = seg_[k].c;
    const double& d = seg_[k].d;

    double h = x-x0;
    return b + 2*c*h + 3*d*h*h;
}

double Spline1D::eval_second_derivative(int k, double x)
{
    // Evaluates the spline segment derivative C_k'' at x
    const double& x0 = seg_[k].x0;
    const double& a = seg_[k].a;
    const double& b = seg_[k].b;
    const double& c = seg_[k].c;
    const double& d = seg_[k].d;

    double h = x-x0;
    return 2*c + 6*d*h;
}