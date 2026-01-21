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

std::vector<double> gauss_elimination(MatrixSystem& sys) {
    int n = sys.n;
    std::vector<double> x(n);

    for (int i = 0; i < n; i++) {
        
        for (int k = i + 1; k < n; k++) {
            double factor = sys.A[k][i] / sys.A[i][i];
            
            for (int j = i; j < n; j++) {
                sys.A[k][j] -= factor * sys.A[i][j];
            }
            sys.r[k] -= factor * sys.r[i];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += sys.A[i][j] * x[j];
        }
        x[i] = (sys.r[i] - sum) / sys.A[i][i];
    }

    return x;
}