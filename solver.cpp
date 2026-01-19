#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cstdlib>

// Problem constants
const double L_val = 2.0;       // Interval length
const double U_0 = 1.0;         // Dirichlet condition u(0) = 1

// Structure storing tridiagonal matrix coefficients
// Matrix equation: A * x = b
struct TridiagonalSystem {
    std::vector<double> a; // lower diagonal
    std::vector<double> b; // main diagonal
    std::vector<double> c; // upper diagonal
    std::vector<double> r; // right-hand side vector (rhs)
    int n;                 // size (number of unknowns)
};

// Shape functions (basis functions) on the reference element [-1, 1] are not explicitly needed here,
// because we integrate over the actual element [x_i, x_{i+1}].
// We use linear shape functions:
// N1 = (x_next - x) / h
// N2 = (x - x_curr) / h

// Gauss-Legendre Quadrature (2 points)
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

int main(int argc, char* argv[]) {
    // 1. Get parameter n
    int n_elements = 0;
    if (argc > 1) {
        n_elements = std::atoi(argv[1]);
    }
    
    if (n_elements < 2) {
        std::cout << "Please provide the number of elements n as a parameter (min. 2). Setting n=10 by default." << std::endl;
        n_elements = 10;
    }

    double h = L_val / n_elements;
    int n_nodes = n_elements + 1;
    
    // Unknowns are w_1, ..., w_n (since w_0 is fixed by Dirichlet = 0 for function w)
    // The system size is n_elements.
    // Vector indexing: 0 corresponds to node x_1, n-1 corresponds to node x_n.
    
    TridiagonalSystem sys;
    sys.n = n_elements;
    sys.a.resize(sys.n, 0.0);
    sys.b.resize(sys.n, 0.0);
    sys.c.resize(sys.n, 0.0);
    sys.r.resize(sys.n, 0.0);

    // 2. Stiffness matrix aggregation (Assembly)
    // Loop over elements
    for (int k = 0; k < n_elements; ++k) {
        double x_left = k * h;
        double x_right = (k + 1) * h;

        // Local node indices in global numbering (0...n)
        int idx_global_1 = k;
        int idx_global_2 = k + 1;

        // Basis functions on the element [x_left, x_right]:
        // N1(x) = (x_right - x) / h   => dN1/dx = -1/h
        // N2(x) = (x - x_left) / h    => dN2/dx =  1/h

        // Local integrals (element stiffness & mass matrices)
        // B(u,v) on element = int (u'v') dx - int (uv) dx
        
        // Calculate integral values analytically or numerically. 
        // Here numerically as per requirements, although it's the same for polynomials.
        
        auto local_stiffness = [&](int i, int j) { // i, j in {1, 2}
            double dNi = (i == 1) ? -1.0/h : 1.0/h;
            double dNj = (j == 1) ? -1.0/h : 1.0/h;
            // integrand for stiffness: u'v'
            auto integrand_stiff = [&](double x) { return dNi * dNj; };
            // integrand for mass: -uv
            auto integrand_mass = [&](double x) {
                 double val_i = (i == 1) ? (x_right - x)/h : (x - x_left)/h;
                 double val_j = (j == 1) ? (x_right - x)/h : (x - x_left)/h;
                 return -1.0 * val_i * val_j;
            };
            return gauss_integrate(x_left, x_right, integrand_stiff) + 
                   gauss_integrate(x_left, x_right, integrand_mass);
        };

        // Load integral L(v) - B(utilde, v)
        // Right hand side = int (sin(x) + 1) * v dx
        auto local_load = [&](int i) {
            auto integrand = [&](double x) {
                double val_v = (i == 1) ? (x_right - x)/h : (x - x_left)/h;
                return (std::sin(x) + 1.0) * val_v;
            };
            return gauss_integrate(x_left, x_right, integrand);
        };

        // Local 2x2 matrix values
        double K11 = local_stiffness(1, 1);
        double K12 = local_stiffness(1, 2);
        double K21 = local_stiffness(2, 1);
        double K22 = local_stiffness(2, 2);

        double F1 = local_load(1);
        double F2 = local_load(2);

        // Adding to global system of equations
        // Remember: we solve for w_1 ... w_n. w_0 is removed (homogeneous Dirichlet).
        // Mapping: Global node `idx` (if > 0) goes to equation `idx - 1`.

        // Node 1 (left node of the element)
        if (idx_global_1 > 0) {
            int row = idx_global_1 - 1;
            sys.b[row] += K11; // Diagonal
            sys.r[row] += F1;
            if (idx_global_2 <= n_elements) { // Right neighbor
                // This is the element above the diagonal for row `row`
                // i.e., sys.c[row]. Corresponds to column row+1.
                sys.c[row] += K12;
            }
        }

        // Node 2 (right node of the element)
        if (idx_global_2 > 0) { // always true for k>=0
            int row = idx_global_2 - 1;
            sys.b[row] += K22; // Diagonal
            sys.r[row] += F2;
            if (idx_global_1 > 0) { // Left neighbor
                // This is the element below the diagonal for row `row`
                // i.e., sys.a[row]. Corresponds to column row-1.
                // In Thomas alg libraries, a[i] is often the element to the left of diagonal in row i.
                // sys.a[row] corresponds to A[row][row-1]
                sys.a[row] += K21;
            }
        }
    }

    // 3. Incorporating boundary conditions (Robin at x=2)
    // Boundary terms from B(u, v): - u(2)v(2) -> in matrix form, subtract 1 from diagonal at last node.
    // Boundary terms from L(v) - B(utilde, v): + 6*v(2) -> add 6 to RHS at last node.
    
    // Last node has index n_elements, corresponding to row n_elements - 1
    int last_row = n_elements - 1;
    sys.b[last_row] += -1.0; 
    sys.r[last_row] += 6.0;

    // 4. Solving the system (Thomas Algorithm)
    // System A * x = r. 
    // a - lower diag, b - main, c - upper.
    // sys.a[0] is unused (nothing to the left of the first element)
    
    std::vector<double> w(n_elements); // Result w_1 ... w_n
    
    // Forward elimination
    std::vector<double> c_prime(n_elements);
    std::vector<double> d_prime(n_elements);

    c_prime[0] = sys.c[0] / sys.b[0];
    d_prime[0] = sys.r[0] / sys.b[0];

    for (int i = 1; i < n_elements; i++) {
        double temp = sys.b[i] - sys.a[i] * c_prime[i-1];
        c_prime[i] = sys.c[i] / temp;
        d_prime[i] = (sys.r[i] - sys.a[i] * d_prime[i-1]) / temp;
    }

    // Back substitution
    w[n_elements - 1] = d_prime[n_elements - 1];
    for (int i = n_elements - 2; i >= 0; i--) {
        w[i] = d_prime[i] - c_prime[i] * w[i+1];
    }

    // 5. Saving results to CSV
    // Remember: u(x) = w(x) + utilde(x) = w(x) + 1.
    // Node 0: u(0) = 1 (w_0 = 0)
    
    std::ofstream file("results.csv");
    file << "x;u(x)" << "\n";
    file << std::fixed << std::setprecision(6);

    // Save node 0
    file << 0.0 << ";" << 1.0 << "\n";

    // Save nodes 1...n
    for (int i = 0; i < n_elements; ++i) {
        double x = (i + 1) * h;
        double u_val = w[i] + U_0; // Restore shift
        file << x << ";" << u_val << "\n";
    }

    file.close();
    std::cout << "Calculations finished for n=" << n_elements << ". Results saved to results.csv" << std::endl;

    return 0;
}