#include "fem_solver.hpp"
#include "math_utils.hpp"
#include <cmath>
#include <vector>

// Internal structure helper
struct TridiagonalSystem {
    std::vector<double> a, b, c, r;
    int n;
};

std::vector<double> solve_fem(int n_elements) {
    double h = L_VAL / n_elements;
    
    // initialize tridiagonal matrix
    TridiagonalSystem sys;
    sys.n = n_elements;
    sys.a.resize(sys.n, 0.0);
    sys.b.resize(sys.n, 0.0);
    sys.c.resize(sys.n, 0.0);
    sys.r.resize(sys.n, 0.0);

    // main loop
    for (int k = 0; k < n_elements; ++k) {
        // calculating element boundaries and global indices
        double x_left = k * h;
        double x_right = (k + 1) * h;
        int idx_global_1 = k;
        int idx_global_2 = k + 1;

        // left side of equation
        auto local_stiffness = [&](int i, int j) {
            double dNi = (i == 1) ? -1.0/h : 1.0/h;
            double dNj = (j == 1) ? -1.0/h : 1.0/h;
            auto integrand_stiff = [&](double x) { return dNi * dNj; };
            auto integrand_mass = [&](double x) {
                 double val_i = (i == 1) ? (x_right - x)/h : (x - x_left)/h;
                 double val_j = (j == 1) ? (x_right - x)/h : (x - x_left)/h;
                 return -1.0 * val_i * val_j;
            };
            return gauss_integrate(x_left, x_right, integrand_stiff) + 
                   gauss_integrate(x_left, x_right, integrand_mass);
        };

        // right side of equation
        auto local_load = [&](int i) {
            auto integrand = [&](double x) {
                double val_v = (i == 1) ? (x_right - x)/h : (x - x_left)/h;
                return (std::sin(x) + 1.0) * val_v;
            };
            return gauss_integrate(x_left, x_right, integrand);
        };

        double K11 = local_stiffness(1, 1); double K12 = local_stiffness(1, 2);
        double K21 = local_stiffness(2, 1); double K22 = local_stiffness(2, 2);
        double F1 = local_load(1); double F2 = local_load(2);


        // adding calcualted values to global system
        if (idx_global_1 > 0) {
            int row = idx_global_1 - 1;
            sys.b[row] += K11; sys.r[row] += F1;
            if (idx_global_2 <= n_elements) sys.c[row] += K12;
        }

        if (idx_global_2 > 0) {
            int row = idx_global_2 - 1;
            sys.b[row] += K22; sys.r[row] += F2;
            if (idx_global_1 > 0) sys.a[row] += K21;
        }
    }

    // Boundary conditions
    int last_row = n_elements - 1;
    sys.b[last_row] += -1.0; 
    sys.r[last_row] += 6.0;

    // Thomas algorithm to solve tridiagonal system
    std::vector<double> w(n_elements);
    std::vector<double> c_prime(n_elements);
    std::vector<double> d_prime(n_elements);

    c_prime[0] = sys.c[0] / sys.b[0];
    d_prime[0] = sys.r[0] / sys.b[0];

    for (int i = 1; i < n_elements; i++) {
        double temp = sys.b[i] - sys.a[i] * c_prime[i-1];
        c_prime[i] = sys.c[i] / temp;
        d_prime[i] = (sys.r[i] - sys.a[i] * d_prime[i-1]) / temp;
    }

    w[n_elements - 1] = d_prime[n_elements - 1];
    for (int i = n_elements - 2; i >= 0; i--) {
        w[i] = d_prime[i] - c_prime[i] * w[i+1];
    }

    return w;
}