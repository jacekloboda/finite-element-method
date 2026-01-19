#pragma once
#include <functional>

// Gauss-Legendre Quadrature (2 points)
double gauss_integrate(double x_start, double x_end, std::function<double(double)> func);