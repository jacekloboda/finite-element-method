#pragma once
#include <vector>

// Problem constants
const double L_VAL = 2.0;       // Interval length
const double U_0 = 1.0;         // Dirichlet condition u(0) = 1

// Solves the differential equation using FEM and returns the vector of results w (excluding w_0)
std::vector<double> solve_fem(int n_elements);