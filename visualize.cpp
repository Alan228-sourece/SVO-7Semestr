#include "inmost.h"
#include <iostream>
#include <functional>
#include <cmath>
#include <chrono>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <numeric>

using namespace INMOST;
using namespace std;

vector<double> inmostSolve(std::function<double(double, double, double, double)> f, 
                          std::function<double(double, double)> dostresh, 
                          double dx, double dy, int n, 
                          double& solve_time, int& iterations);

double u(double x, double y) {
    return sin(5.0 * x) * sin(5.0 * y);
}

double f(double x, double y, double dx, double dy) {
    return 25 * (dx + dy) * sin(5.0 * x) * sin(5.0 * y);
}

int main() {
    vector<int> grids = {10, 20, 40, 80, 160, 320};
    vector<double> h_values;
    vector<double> maxers;
    vector<double> l2errs;
    vector<double> solve_times;
    vector<int> iterations_list;
    
    vector<vector<double>> Z_fine;
    int finest_n = 40; 
    
    
    for (size_t idx = 0; idx < grids.size(); idx++) {
        int n = grids[idx];
        double h = 1.0 / n;
        
        cout << "\nСетка " << n << "x" << n << " (h = " << h << ")" << endl;
        
        // Решаем задачу
        double solve_time = 0.0;
        int iterations = 0;
        auto start_time = chrono::high_resolution_clock::now();
        
        vector<double> res = inmostSolve(f, u, 1.0, 1.0, n, solve_time, iterations);
        
        auto end_time = chrono::high_resolution_clock::now();
        double total_time = chrono::duration<double>(end_time - start_time).count();
        
        vector<vector<double>> Z(n + 1, vector<double>(n + 1, 0));
        
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                int idx_res = (i - 1) * (n - 1) + (j - 1);
                Z[i][j] = res[idx_res];
            }
        }
        
        // Граничные условия
        for (int i = 0; i < n + 1; i++) {
            Z[0][i] = u(0.0, i * h);
            Z[n][i] = u(1.0, i * h);
            Z[i][0] = u(i * h, 0.0);
            Z[i][n] = u(i * h, 1.0);
        }
        
        if (n == finest_n) {
            Z_fine = Z;
        }
        
        double max_error = 0.0;
        double l2_error = 0.0;
        int interior_points = (n - 1) * (n - 1);
        
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                double exact = u(i * h, j * h);
                double error = fabs(Z[i][j] - exact);
                max_error = max(max_error, error);
                l2_error += error * error;
            }
        }
        
        l2_error = sqrt(l2_error / interior_points);
        
        h_values.push_back(h);
        maxers.push_back(max_error);
        l2errs.push_back(l2_error);
        solve_times.push_back(total_time);
        iterations_list.push_back(iterations);
    }
    
    if (!Z_fine.empty()) {
        ofstream grid_file("solution_grid.txt");
        if (grid_file.is_open()) {
            for (size_t i = 0; i < Z_fine.size(); i++) {
                for (size_t j = 0; j < Z_fine[i].size(); j++) {
                    grid_file << Z_fine[i][j] << " ";
                }
                grid_file << endl;
            }
            grid_file.close();
            cout << "\nДанные решения сохранены в solution_grid.txt" << endl;
        }
    }

    ofstream conv_file("convergence_data.txt");
    if (conv_file.is_open()) {
        conv_file << "# h max_error l2_error solve_time iterations" << endl;
        for (size_t i = 0; i < grids.size(); i++) {
            conv_file << scientific << setprecision(6) 
                     << h_values[i] << " "
                     << maxers[i] << " "
                     << l2errs[i] << " "
                     << solve_times[i] << " "
                     << iterations_list[i] << endl;
        }
        conv_file.close();
    }
    
    ofstream py_script("visualize_all.py");
    if (py_script.is_open()) {
        py_script << R"(
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
try:
    n = 40
    with open('solution_grid.txt', 'r') as f:
        lines = f.readlines()

    Z = np.zeros((n+1, n+1))
    for i in range(n+1):
        Z[i, :] = list(map(float, lines[i].strip().split()))

    x = np.linspace(0, 1, n+1)
    y = np.linspace(0, 1, n+1)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(14, 6))
    
    # 3D поверхность
    ax1 = fig.add_subplot(121, projection='3d')
    surf = ax1.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('u(x,y)')
    ax1.set_title('Численное решение уравнения диффузии')
    fig.colorbar(surf, ax=ax1, shrink=0.5)
    
    # Контурный график
    ax2 = fig.add_subplot(122)
    contour = ax2.contourf(X, Y, Z, 20, cmap='viridis')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_title('Контуры решения (n=40)')
    fig.colorbar(contour, ax=ax2, shrink=0.5)
    
    plt.tight_layout()
    plt.savefig('solution_plot.png', dpi=150)
    plt.show()
    print("  График решения сохранен как solution_plot.png")
except Exception as e:
    print(f"  Ошибка при построении решения: {e}")

print("\n2. Построение графиков сходимости...")
try:
    # Чтение данных
    data = np.loadtxt('convergence_data.txt', comments='#')
    h = data[:, 0]
    max_err = data[:, 1]
    l2_err = data[:, 2]
    times = data[:, 3]
    iterations = data[:, 4]
    
    # График 1: Зависимость ошибок от шага сетки
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Линейно-логарифмический масштаб
    ax1.loglog(h, max_err, 'o-', linewidth=2, markersize=8, label='Max ошибка')
    ax1.loglog(h, l2_err, 's-', linewidth=2, markersize=8, label='L2 ошибка')
    
    # Добавляем линии сходимости
    order1 = max_err[0] * (h / h[0])**1
    order2 = max_err[0] * (h / h[0])**2
    ax1.loglog(h, order1, 'k--', linewidth=1, label='O(h)')
    ax1.loglog(h, order2, 'k:', linewidth=1, label='O(h²)')
    
    ax1.set_xlabel('Шаг сетки h', fontsize=12)
    ax1.set_ylabel('Норма ошибки', fontsize=12)
    ax1.set_title('Зависимость ошибок от шага сетки', fontsize=14)
    ax1.grid(True, which="both", ls="-", alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Обычный масштаб
    ax2.plot(h, max_err, 'o-', linewidth=2, markersize=8, label='Max ошибка')
    ax2.plot(h, l2_err, 's-', linewidth=2, markersize=8, label='L2 ошибка')
    ax2.set_xlabel('Шаг сетки h', fontsize=12)
    ax2.set_ylabel('Норма ошибки', fontsize=12)
    ax2.set_title('Зависимость ошибок от шага сетки', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('convergence_error.png', dpi=150)
    plt.show()
    print("  График ошибок сохранен как convergence_error.png")
    
    # График 2: Зависимость времени решения от размера сетки
    n_points = 1.0 / h
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    ax1.plot(n_points, times, 'o-', linewidth=2, markersize=8, color='red')
    ax1.set_xlabel('Число узлов (1/h)', fontsize=12)
    ax1.set_ylabel('Время решения, с', fontsize=12)
    ax1.set_title('Зависимость времени решения от размера сетки', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Логарифмический масштаб для определения порядка
    ax2.loglog(n_points, times, 'o-', linewidth=2, markersize=8, color='red')
    ax2.set_xlabel('Число узлов (1/h)', fontsize=12)
    ax2.set_ylabel('Время решения, с', fontsize=12)
    ax2.set_title('Зависимость времени решения (лог. масштаб)', fontsize=14)
    ax2.grid(True, which="both", ls="-", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('convergence_time.png', dpi=150)
    plt.show()
    print("  График времени сохранен как convergence_time.png")
    
    # График 3: Зависимость числа итераций от размера сетки
    fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    ax1.plot(n_points, iterations, 'o-', linewidth=2, markersize=8, color='green')
    ax1.set_xlabel('Число узлов (1/h)', fontsize=12)
    ax1.set_ylabel('Число итераций', fontsize=12)
    ax1.set_title('Зависимость числа итераций от размера сетки', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    ax2.loglog(n_points, iterations, 'o-', linewidth=2, markersize=8, color='green')
    ax2.set_xlabel('Число узлов (1/h)', fontsize=12)
    ax2.set_ylabel('Число итераций', fontsize=12)
    ax2.set_title('Зависимость числа итераций (лог. масштаб)', fontsize=14)
    ax2.grid(True, which="both", ls="-", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('convergence_iterations.png', dpi=150)
    plt.show()
    print("  График итераций сохранен как convergence_iterations.png")
    
    # Сводная таблица сходимости
    print("\n3. Таблица сходимости:")
    print("┌──────┬──────────┬──────────────┬──────────────┬───────────┬────────────┐")
    print("│   N  │     h    │  Max ошибка  │   L2 ошибка  │  Время,с  │  Итераций  │")
    print("├──────┼──────────┼──────────────┼──────────────┼───────────┼────────────┤")
    for i in range(len(h)):
        n_size = int(1.0/h[i])
        print(f"│ {n_size:4d} │ {h[i]:8.6f} │ {max_err[i]:12.2e} │ {l2_err[i]:12.2e} │ {times[i]:9.4f} │ {int(iterations[i]):10d} │")
    print("└──────┴──────────┴──────────────┴──────────────┴───────────┴────────────┘")
    
    # Расчет порядка сходимости
    print("\n4. Порядок сходимости:")
    for i in range(1, len(h)):
        if max_err[i-1] > 0 and max_err[i] > 0:
            p_max = np.log(max_err[i-1]/max_err[i]) / np.log(h[i-1]/h[i])
            p_l2 = np.log(l2_err[i-1]/l2_err[i]) / np.log(h[i-1]/h[i])
            print(f"  h={h[i-1]:.4f} -> h={h[i]:.4f}: p_max = {p_max:.3f}, p_l2 = {p_l2:.3f}")
    
except Exception as e:
    print(f"  Ошибка при построении графиков сходимости: {e}")

print("\nВсе графики построены!")
print("Файлы:")
print("  - solution_plot.png - визуализация решения")
print("  - convergence_error.png - график ошибок")
print("  - convergence_time.png - график времени")
print("  - convergence_iterations.png - график итераций")
)";
        py_script.close();
    }
    
        system("python visualize_all.py");
    
    return 0;
}