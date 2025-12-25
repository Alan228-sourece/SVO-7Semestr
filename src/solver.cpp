#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <functional>
#include "inmost.h"

using namespace INMOST;
using namespace std;

void fillM(Sparse::Matrix& A, Sparse::Vector& b, 
           double dx, double dy, int n,
           std::function<double(double, double, double, double)> f,
           std::function<double(double, double)> dostresh) {
    
    double h = 1.0 / n;
    int inpnt = (n - 1) * (n - 1);
    
    auto index = [n](int i, int j) {
        return (i-1) * (n-1) + (j-1);
    };
    
    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < n; ++j) {
            int idx = index(i, j);
            double xi = i * h;
            double yj = j * h;
            A[idx][idx] = 2.0 * (dx + dy) / (h * h);
            if (i > 1) {
                int lidx = index(i-1, j);
                A[idx][lidx] = -dx / (h * h);
            }
            if (i < n-1) {
                int ridx = index(i+1, j);
                A[idx][ridx] = -dx / (h * h);
            }
            if (j > 1) {
                int bidx = index(i, j-1);
                A[idx][bidx] = -dy / (h * h);
            }
            if (j < n-1) {
                int tidx = index(i, j+1);
                A[idx][tidx] = -dy / (h * h);
            }
            b[idx] = f(xi, yj, dx, dy);
            
            // граничные условия
            if (i == 1) {
                b[idx] += (dx / (h * h)) * dostresh(0.0, yj);
            }
            if (i == n-1) {
                b[idx] += (dx / (h * h)) * dostresh(1.0, yj);
            }
            if (j == 1) {
                b[idx] += (dy / (h * h)) * dostresh(xi, 0.0);
            }
            if (j == n-1) {
                b[idx] += (dy / (h * h)) * dostresh(xi, 1.0);
            }
        }
    }
}

vector<double> inmostSolve(std::function<double(double, double, double, double)> f, 
                          std::function<double(double, double)> dostresh, 
                          double dx, double dy, int n,
                          double& soltime, int& iterations) {
    
    int inpnt = (n - 1) * (n - 1);
    vector<double> res(inpnt);
    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector x;
    A.SetInterval(0, inpnt);
    b.SetInterval(0, inpnt);
    x.SetInterval(0, inpnt);
    fillM(A, b, dx, dy, n, f, dostresh);
    Solver solver("inner_mptilu2");
    solver.SetParameter("absolute_tolerance", "1e-10");
    solver.SetParameter("relative_tolerance", "1e-10");
    solver.SetParameter("drop_tolerance", "0.001");
    
    solver.SetMatrix(A);
     solver.Solve(b, x);
    iterations = solver.Iterations();
    
   
    for (int i = 0; i < inpnt; ++i) {
        res[i] = x[i];
    }
    
    return res;
}