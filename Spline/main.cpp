#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "CubicSpline.hpp"

/** Функция для вычисления истинного значения y = e^x **/
double true_function(double x) {
    return std::exp(x);
}

/** Функция для вычисления второй производной y'' = e^x **/
double true_second_derivative(double x) {
    return std::exp(x);
}

int main() {
    std::vector<int> N_values = { 5, 10, 20, 40, 80, 160 };
    std::vector<double> errors_natural;
    std::vector<double> errors_true;

    // Открытие файлов для записи ошибок
    std::ofstream natural_error_file("natural_spline_error.txt");
    std::ofstream true_error_file("true_spline_error.txt");

    for (auto N : N_values) {
        // Генерация узлов
        std::vector<double> x(N);
        std::vector<double> y(N);
        double a = 0.0;
        double b = 10.0;
        double step = (b - a) / (N - 1);
        for (int i = 0; i < N; ++i) {
            x[i] = a + i * step;
            y[i] = true_function(x[i]);
        }

        // Естественный сплайн (границы вторых производных равны 0)
        CubicSpline<double, double> natural_spline(x, y, 0.0, 0.0);

        // Сплайн с истинными вторыми производными на границах
        double fpp_a = true_second_derivative(a);
        double fpp_b = true_second_derivative(b);
        CubicSpline<double, double> true_spline(x, y, fpp_a, fpp_b);

        // Оценка ошибки
        int M = 1000;
        double error_natural = 0.0;
        double error_true = 0.0;
        double step_eval = (b - a) / (M - 1);
        for (int i = 0; i < M; ++i) {
            double xi = a + i * step_eval;
            double yi_true = true_function(xi);
            double yi_natural = natural_spline.interpolate(xi);
            double yi_true_spline = true_spline.interpolate(xi);
            double current_error_natural = std::abs(yi_true - yi_natural);
            double current_error_true = std::abs(yi_true - yi_true_spline);
            if (current_error_natural > error_natural) {
                error_natural = current_error_natural;
            }
            if (current_error_true > error_true) {
                error_true = current_error_true;
            }
        }

        errors_natural.push_back(error_natural);
        errors_true.push_back(error_true);

        // Запись в файлы
        natural_error_file << std::log(N) << " " << std::log(error_natural) << "\n";
        true_error_file << std::log(N) << " " << std::log(error_true) << "\n";

        std::cout << "N = " << N
            << ", Error Natural = " << error_natural
            << ", Error True = " << error_true << std::endl;
    }

    natural_error_file.close();
    true_error_file.close();

    // Здесь можно использовать внешние инструменты для построения графиков, например, Python с matplotlib

    return 0;
}
