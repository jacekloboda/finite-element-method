#include "math_utils.hpp"
#include <cmath>

double gauss_integrate(double x_start, double x_end, std::function<double(double)> func) {
    // Points and weights for the interval [-1, 1]
    const double xi1 = -1.0 / std::sqrt(3.0);
    const double xi2 =  1.0 / std::sqrt(3.0);
    const double w = 1.0;

    // Transformation to [x_start, x_end]
    double center = (x_start + x_end) / 2.0;
    double half_width = (x_end - x_start) / 2.0;

    double val1 = func(center + half_width * xi1);
    double val2 = func(center + half_width * xi2);

    return half_width * (w * val1 + w * val2);
}