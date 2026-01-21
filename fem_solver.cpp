#include "fem_solver.hpp"
#include "math_utils.hpp"
#include <cmath>
#include <vector>
#include <iostream>


std::vector<double> solve_fem(int n_elements) {
    double h = L_VAL / n_elements;
    
    MatrixSystem sys;
    sys.n = n_elements;
    
    sys.A.resize(sys.n, std::vector<double>(sys.n, 0.0));
    sys.r.resize(sys.n, 0.0);

    // main loop
    for (int k = 0; k < n_elements; ++k) {
        // Element boundaries and global node indices
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
        

        if (idx_global_1 > 0) {
            int row = idx_global_1 - 1;
            
            sys.A[row][row] += K11;
            sys.r[row]      += F1;

            if (idx_global_2 <= n_elements) {
                int col = idx_global_2 - 1;
                sys.A[row][col] += K12;
            }
        }


        if (idx_global_2 > 0) {
            int row = idx_global_2 - 1;
            
            sys.A[row][row] += K22;
            sys.r[row]      += F2;

            if (idx_global_1 > 0) {
                int col = idx_global_1 - 1;
                sys.A[row][col] += K21;
            }
        }
    }

    int last_row = n_elements - 1;
    sys.A[last_row][last_row] += -1.0; 
    sys.r[last_row] += 6.0;

    return gauss_elimination(sys);
}