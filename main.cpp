#include "inmost.h"
#include <iostream>
#include <functional>
#include <cmath>
#include <chrono>
#include <vector>

using namespace INMOST;
using namespace std;

vector<double> inmostSolve(std::function<double(double, double, double, double)> f, 
                          std::function<double(double, double)> dostresh, 
                          double dx, double dy, int n, 
                          double& solve_time, int& iterations);
double u(double x, double y) {
    return sin(5.0 * x) * cos(5.0 * y);
}

double f(double x, double y, double dx, double dy) {
    return 25.0 * (dx + dy) * sin(5.0 * x) * cos(5.0 * y);
}

void fillM(Sparse::Matrix& A, Sparse::Vector& b, 
           double dx, double dy,
           int n,
           std::function<double(double, double, double, double)> f);


double norm_L2_01(const vector<double>& x, std::function<double(double, double)> u, size_t n) {
    double h = 1.0 / n;
    double norm = 0.0;
    for (size_t i = 1; i < n; i++) {
        for (size_t j = 1; j < n; j++) {
            size_t idx = (i - 1) * (n - 1) + j - 1;
            double diff = u(i * h, j * h) - x[idx];
            norm += diff * diff;
        }
    }
    return sqrt(norm) * h;
}

double norm_C_01(const vector<double>& x, std::function<double(double, double)> u, size_t n) {
    double h = 1.0 / n;
    double norm = 0.0;
    for (size_t i = 1; i < n; i++) {
        for (size_t j = 1; j < n; j++) {
            size_t idx = (i - 1) * (n - 1) + j - 1;
            norm = std::max(norm, abs(u(i * h, j * h) - x[idx]));
        }
    }
    return norm;
}

int main(int argc, char *argv[]) {
  /*    setlocale(LC_ALL, "");

    size_t max_n = 600;
    cout<<"NIGGER";
    for (size_t n = 2; n < max_n; n *= 2) {
        std::cout << n << ' ';
        auto start = std::chrono::steady_clock::now();
        vector<double> x = inmostSolve(f,u, 1.0, 1.0, n);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> time = end - start;
        std::cout << time.count() << ' ' << norm_L2_01(x, u, n) << ' ' << norm_C_01(x, u, n) << std::endl;
    }
    return 0;*/
}