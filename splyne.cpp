#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <type_traits>

// Определение класса для работы с трехдиагональной матрицей
template<typename Type>
class ThreeDiagonalMatrix {
public:
    std::vector<Type> a; // Нижняя диагональ (a_1, a_2, ..., a_{n-1})
    std::vector<Type> b; // Главная диагональ (b_0, b_1, ..., b_{n-1})
    std::vector<Type> c; // Верхняя диагональ (c_0, c_1, ..., c_{n-2})

    // Конструктор
    ThreeDiagonalMatrix() {}

    // Инициализация матрицы с размером n
    ThreeDiagonalMatrix(size_t n) {
        if (n == 0) throw std::invalid_argument("Размер матрицы должен быть больше нуля.");
        a.resize(n - 1, 0);
        b.resize(n, 0);
        c.resize(n - 1, 0);
    }

    // Метод для вывода матрицы (для отладки)
    void print() const {
        size_t n = b.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (j == i - 1) {
                    std::cout << a[i - 1] << " ";
                }
                else if (j == i) {
                    std::cout << b[i] << " ";
                }
                else if (j == i + 1) {
                    std::cout << c[i] << " ";
                }
                else {
                    std::cout << "0 ";
                }
            }
            std::cout << std::endl;
        }
    }
};

// Определение типов для деления и разности
template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

// Функция для решения системы уравнений методом прогонки (метод Томаса)
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType>& matrix,
                                           const std::vector<cType>& column) {
    size_t n = matrix.b.size();
    if (column.size() != n) {
        throw std::invalid_argument("Размер правой части должен соответствовать размеру матрицы.");
    }

    // Векторы для модифицированных коэффициентов
    std::vector<DivisType<cType, mType>> c_star(n - 1);
    std::vector<DivisType<cType, mType>> d_star(n);

    // Прямой ход
    c_star[0] = matrix.c[0] / matrix.b[0];
    d_star[0] = column[0] / matrix.b[0];
    for (size_t i = 1; i < n - 1; ++i) {
        DivisType<cType, mType> m = matrix.b[i] - matrix.a[i - 1] * c_star[i - 1];
        c_star[i] = matrix.c[i] / m;
        d_star[i] = (column[i] - matrix.a[i - 1] * d_star[i - 1]) / m;
    }
    // Последний элемент d_star
    d_star[n - 1] = (column[n - 1] - matrix.a[n - 2] * d_star[n - 2]) / (matrix.b[n - 1] - matrix.a[n - 2] * c_star[n - 2]);

    // Обратный ход
    std::vector<DivisType<cType, mType>> x(n);
    x[n - 1] = d_star[n - 1];
    for (size_t i = n - 1; i-- > 0;) {
        x[i] = d_star[i] - matrix.c[i] * x[i + 1];
    }

    return x;
}

// Класс кубического сплайна
template<typename xType, typename yType>
class CubicSpline {
private:
    std::vector<xType> x;          // Узлы интерполяции
    std::vector<yType> y;          // Значения функции в узлах
    std::vector<yType> M;          // Вторые производные сплайна
    bool natural;                   // Флаг для естественных сплайнов

public:
    /**
     * Конструктор кубического сплайна.
     *
     * @param points Вектор значений x.
     * @param values Вектор значений y.
     * @param first  Значение второй производной в левом конце (для сплайнов с заданными производными).
     *               Для естественных сплайнов установить значение natural=true.
     * @param second Значение второй производной в правом конце (для сплайнов с заданными производными).
     *               Для естественных сплайнов установить значение natural=true.
     * @param isNatural Флаг, указывающий на тип граничных условий.
     */
    CubicSpline(const std::vector<xType>& points,
               const std::vector<yType>& values,
               const yType& first = 0,
               const yType& second = 0,
               bool isNatural = true)
        : x(points), y(values), natural(isNatural) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("Размеры векторов x и y должны совпадать.");
        }
        if (x.size() < 2) {
            throw std::invalid_argument("Должно быть как минимум две точки для интерполяции.");
        }

        computeSecondDerivatives(first, second);
    }

    /**
     * Метод для вычисления вторых производных сплайна.
     *
     * @param first  Значение второй производной в левом конце.
     * @param second Значение второй производной в правом конце.
     */
    void computeSecondDerivatives(const yType& first, const yType& second) {
        size_t n = x.size();
        std::vector<yType> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            h[i] = x[i + 1] - x[i];
            if (h[i] <= 0) {
                throw std::invalid_argument("Значения x должны быть строго возрастающими.");
            }
        }

        // Формирование системы для метода прогонки
        ThreeDiagonalMatrix<yType> matrix(n);
        std::vector<yType> d(n, 0);

        if (natural) {
            // Естественные сплайны: M0 = Mn = 0
            matrix.b[0] = 1;
            d[0] = 0;

            for (size_t i = 1; i < n - 1; ++i) {
                matrix.a[i - 1] = h[i - 1];
                matrix.b[i] = 2 * (h[i - 1] + h[i]);
                matrix.c[i] = h[i];
                d[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
            }

            matrix.b[n - 1] = 1;
            d[n - 1] = 0;
        }
        else {
            // Сплайны с заданными вторыми производными на концах: M0 = first, Mn = second
            matrix.b[0] = 1;
            d[0] = first;

            for (size_t i = 1; i < n - 1; ++i) {
                matrix.a[i - 1] = h[i - 1];
                matrix.b[i] = 2 * (h[i - 1] + h[i]);
                matrix.c[i] = h[i];
                d[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
            }

            matrix.b[n - 1] = 1;
            d[n - 1] = second;
        }

        // Решение системы уравнений для нахождения вторых производных
        M = solve(matrix, d);
    }

    /**
     * Метод для интерполяции значения функции в заданной точке.
     *
     * @param xi Точка, в которой необходимо интерполировать.
     * @return Интерполированное значение.
     */
    yType interpolate(const xType& xi) const noexcept {
        size_t n = x.size();

        // Поиск интервала, содержащего xi
        size_t i = 0;
        if (xi <= x[0]) {
            i = 0;
        }
        else if (xi >= x[n - 1]) {
            i = n - 2;
        }
        else {
            // Используем бинарный поиск для эффективности
            size_t low = 0;
            size_t high = n - 1;
            while (low <= high) {
                size_t mid = low + (high - low) / 2;
                if (x[mid] <= xi && xi < x[mid + 1]) {
                    i = mid;
                    break;
                }
                else if (xi < x[mid]) {
                    if (mid == 0) break;
                    high = mid - 1;
                }
                else {
                    low = mid + 1;
                }
            }
        }

        yType hi = x[i + 1] - x[i];
        yType a = (x[i + 1] - xi) / hi;
        yType b = (xi - x[i]) / hi;
        yType ai = y[i];
        yType bi = y[i + 1];
        yType Mi = M[i];
        yType Mi1 = M[i + 1];

        // Кубическая формула сплайна
        yType S = a * ai + b * bi + ((a * a * a - a) * Mi + (b * b * b - b) * Mi1) * (hi * hi) / 6;
        return S;
    }
};

// Функция для генерации равномерных узлов на отрезке [a, b]
std::vector<double> generateUniformNodes(double a, double b, size_t N) {
    std::vector<double> nodes(N);
    double h = (b - a) / (N - 1);
    for (size_t i = 0; i < N; ++i) {
        nodes[i] = a + i * h;
    }
    return nodes;
}

int main() {
    // Настройки
    std::vector<size_t> N_values = {5, 10, 20, 40, 80, 160};
    double a = 0.0;
    double b = 10.0;
    size_t num_test_points = 1000; // Количество точек для оценки ошибки

    // Файлы для сохранения результатов
    std::ofstream naturalFile("natural_spline_errors.csv");
    std::ofstream exactFile("exact_spline_errors.csv");

    if (!naturalFile.is_open() || !exactFile.is_open()) {
        std::cerr << "Не удалось открыть файл для записи результатов." << std::endl;
        return 1;
    }

    // Запись заголовков
    naturalFile << "N,MaxError\n";
    exactFile << "N,MaxError\n";

    naturalFile << std::fixed << std::setprecision(10);
    exactFile << std::fixed << std::setprecision(10);

    // Функция для генерации y = e^x
    auto func = [](double x) -> double {
        return std::exp(x);
    };

    // Функция для генерации второй производной y'' = e^x
    auto second_deriv = [](double x) -> double {
        return std::exp(x);
    };

    // Цикл по различным значениям N
    for (size_t N : N_values) {
        // Генерация узлов
        std::vector<double> nodes = generateUniformNodes(a, b, N);

        // Вычисление значений функции в узлах
        std::vector<double> values(N);
        for (size_t i = 0; i < N; ++i) {
            values[i] = func(nodes[i]);
        }

        // === Естественные сплайны ===
        // Граничные условия: M0 = Mn = 0
        CubicSpline<double, double> naturalSpline(nodes, values, 0.0, 0.0, true);

        // Оценка ошибки
        double maxErrorNatural = 0.0;
        double step = (b - a) / (num_test_points - 1);
        for (size_t i = 0; i < num_test_points; ++i) {
            double xi = a + i * step;
            double yi_true = func(xi);
            double yi_spline = naturalSpline.interpolate(xi);
            double error = std::abs(yi_true - yi_spline);
            if (error > maxErrorNatural) {
                maxErrorNatural = error;
            }
        }

        // Запись результата
        naturalFile << N << "," << maxErrorNatural << "\n";

        // === Сплайны с заданными вторыми производными ===
        // Граничные условия: M0 = y''(a), Mn = y''(b)
        double M0 = second_deriv(a);
        double Mn = second_deriv(b);
        CubicSpline<double, double> exactSpline(nodes, values, M0, Mn, false);

        // Оценка ошибки
        double maxErrorExact = 0.0;
        for (size_t i = 0; i < num_test_points; ++i) {
            double xi = a + i * step;
            double yi_true = func(xi);
            double yi_spline = exactSpline.interpolate(xi);
            double error = std::abs(yi_true - yi_spline);
            if (error > maxErrorExact) {
                maxErrorExact = error;
            }
        }

        // Запись результата
        exactFile << N << "," << maxErrorExact << "\n";

        // Вывод в консоль для отслеживания процесса
        std::cout << "N = " << N << " | Natural Spline Max Error = " << maxErrorNatural
                  << " | Exact Spline Max Error = " << maxErrorExact << std::endl;
    }

    // Закрытие файлов
    naturalFile.close();
    exactFile.close();

    std::cout << "Расчеты завершены. Результаты сохранены в файлы natural_spline_errors.csv и exact_spline_errors.csv." << std::endl;

    return 0;
}
