#pragma once
#include <functional>

struct MatrixSystem {
    std::vector<std::vector<double>> A;
    std::vector<double> r;              
    int n;
};

double gauss_integrate(double x_start, double x_end, std::function<double(double)> func);

std::vector<double> gauss_elimination(MatrixSystem& sys);