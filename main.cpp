#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>

#include "fem_solver.hpp"

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

    // 2. Solve the system
    std::vector<double> w = solve_fem(n_elements);

    // 3. Save results to CSV
    std::ofstream file("results.csv");
    file << "x,u(x)" << "\n";
    file << std::fixed << std::setprecision(6);

    // Save node 0
    file << 0.0 << "," << U_0 << "\n";

    double h = L_VAL / n_elements;
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