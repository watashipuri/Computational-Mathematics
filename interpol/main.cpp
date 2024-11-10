#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include "header.h"

int main() {
    // Функция для интерполяции
    auto func = [](double x) { return std::exp(x); };
    auto deriv = [](double x) { return std::exp(x); };

    // Массивы для количества узлов и длин отрезков
    std::vector<unsigned int> N_values = { 3, 4, 5 };
    std::vector<double> L_values = { 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625 };

    // Открытие файла для записи результатов
    std::ofstream outFile("results.csv");
    outFile << "N,L,Max Error\n";  // Заголовки для CSV

    // Основной цикл по значениям N и L
    for (unsigned int N : N_values) {
        for (double L : L_values) {
            // Создаём массивы для узлов, значений функции и производных
            std::array<double, 5> points; // Максимальное N=5
            std::array<double, 5> values;
            std::array<double, 5> derivatives;

            // Заполняем массивы для текущего значения N
            for (unsigned int i = 0; i < N; ++i) {
                points[i] = L * i / (N - 1);  // Равномерное распределение узлов
                values[i] = func(points[i]);
                derivatives[i] = deriv(points[i]);
            }

            // Переменные для хранения максимальной ошибки
            double max_error = 0.0;

            // Создание объекта интерполятора на основе значения N
            if (N == 3) {
                NewtonInterpolator<double, double, 3> interpolator({ points[0], points[1], points[2] },
                    { values[0], values[1], values[2] },
                    { derivatives[0], derivatives[1], derivatives[2] });
                // Вычисляем ошибку интерполяции в 1000 точках
                for (int i = 0; i <= 1000; ++i) {
                    double x = L * i / 1000.0;
                    double interp_value = interpolator.interpolate(x);
                    double exact_value = func(x);
                    double error = std::abs(interp_value - exact_value);
                    if (error > max_error) {
                        max_error = error;
                    }
                }
            }
            else if (N == 4) {
                NewtonInterpolator<double, double, 4> interpolator({ points[0], points[1], points[2], points[3] },
                    { values[0], values[1], values[2], values[3] },
                    { derivatives[0], derivatives[1], derivatives[2], derivatives[3] });
                // Вычисляем ошибку интерполяции в 1000 точках
                for (int i = 0; i <= 1000; ++i) {
                    double x = L * i / 1000.0;
                    double interp_value = interpolator.interpolate(x);
                    double exact_value = func(x);
                    double error = std::abs(interp_value - exact_value);
                    if (error > max_error) {
                        max_error = error;
                    }
                }
            }
            else if (N == 5) {
                NewtonInterpolator<double, double, 5> interpolator({ points[0], points[1], points[2], points[3], points[4] },
                    { values[0], values[1], values[2], values[3], values[4] },
                    { derivatives[0], derivatives[1], derivatives[2], derivatives[3], derivatives[4] });
                // Вычисляем ошибку интерполяции в 1000 точках
                for (int i = 0; i <= 1000; ++i) {
                    double x = L * i / 1000.0;
                    double interp_value = interpolator.interpolate(x);
                    double exact_value = func(x);
                    double error = std::abs(interp_value - exact_value);
                    if (error > max_error) {
                        max_error = error;
                    }
                }
            }

            // Запись результатов для текущих N и L в файл
            outFile << N << "," << L << "," << max_error << "\n";
        }
    }

    outFile.close();
    std::cout << "Results saved as results.csv" << std::endl;
    return 0;
}
