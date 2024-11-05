#pragma once

#include <iostream>
#include <array>
#include <eigen-3.4.0\Eigen\Dense>

template<typename T, unsigned int N>
struct DiffCoefficients {
    T mainCoefficient;
    std::array<T, N> additionalCoefficients;
};

template<typename T, unsigned int N>
DiffCoefficients<T, N> computeDiffCoefficients(const std::array<T, N>& nodes) noexcept {
    Eigen::Matrix<T, N + 1, N + 1> matrixA;

    for (unsigned int col = 0; col <= N; ++col) {
        matrixA(0, col) = static_cast<T>(1);
    }

    for (unsigned int row = 1; row <= N; ++row) {
        matrixA(row, 0) = static_cast<T>(0);
        for (unsigned int col = 1; col <= N; ++col) {
            matrixA(row, col) = matrixA(row - 1, col) * nodes[col - 1];
        }
    }

    Eigen::Matrix<T, N + 1, 1> vectorB = Eigen::Matrix<T, N + 1, 1>::Zero();
    vectorB(1) = static_cast<T>(1);

    Eigen::Matrix<T, N + 1, 1> solution = matrixA.colPivHouseholderQr().solve(vectorB);

    T mainCoef = solution(0);

    std::array<T, N> otherCoefs;
    for (unsigned int i = 0; i < N; ++i) {
        otherCoefs[i] = solution(i + 1);
    }

    return DiffCoefficients<T, N>{ mainCoef, otherCoefs };
}

template<typename T, unsigned int N>
T approximateDerivative(const T x0, const T h, const std::array<T, N>& nodes) {
    DiffCoefficients<T, N> coeffs = computeDiffCoefficients<T, N>(nodes);

    T result = coeffs.mainCoefficient * std::exp(x0) / h;
    for (unsigned int i = 0; i < N; ++i) {
        result += coeffs.additionalCoefficients[i] * std::exp(x0 + nodes[i] * h) / h;
    }

    return result;
}
