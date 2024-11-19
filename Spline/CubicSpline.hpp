#ifndef CUBICSPLINE_HPP
#define CUBICSPLINE_HPP

#include <vector>
#include <type_traits>
#include <cmath>
#include <stdexcept>

/** Класс для работы с трехдиагональной матрицей */
template<typename Type>
class ThreeDiagonalMatrix {
public:
    std::vector<Type> lower; // Нижняя диагональ (a)
    std::vector<Type> main;  // Главная диагональ (b)
    std::vector<Type> upper; // Верхняя диагональ (c)

    ThreeDiagonalMatrix(int size) : lower(size - 1, 0), main(size, 0), upper(size - 1, 0) {}
};

/** Определение типов для деления и разности */
template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

/** Функция для решения методом прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType>& matrix,
    const std::vector<cType>& column);

/**
 * xType - тип аргумента x.
 * yType - тип значения функции y
 */
template<typename xType, typename yType>
class CubicSpline {
    using DeltaXType = DiffType<xType>;
    using DerivType = DivisType<DiffType<yType>, DeltaXType>;
    using Deriv2Type = DivisType<DiffType<DerivType>, DeltaXType>;

    std::vector<xType> points;
    std::vector<yType> values;
    std::vector<Deriv2Type> secondDerivatives;

public:
    CubicSpline(const std::vector<xType>& points, // Значения x
        const std::vector<yType>& values, // Значения y
        const Deriv2Type& first,          // Значение для левой второй производной
        const Deriv2Type& second);        // Значение для правой второй производной

    yType interpolate(const xType& x) const noexcept;
};

#include "CubicSpline.tpp" // Включение реализации шаблонов

#endif // CUBICSPLINE_HPP
