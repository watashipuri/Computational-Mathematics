#pragma once
#include <array>

template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
    std::array<xType, 2 * N> duplicated_x;
    std::array<yType, 2 * N> duplicated_f;
    std::array<std::array<yType, 2 * N>, 2 * N> divided_diff;
    std::array<yType, 2 * N> coefficients;

public:
    NewtonInterpolator(const std::array<xType, N>& points, const std::array<yType, N>& values, const std::array<yType, N>& deriv) noexcept;

    yType interpolate(const xType& x) const noexcept;
};

template<typename xType, typename yType, unsigned int N>
NewtonInterpolator<xType, yType, N>::NewtonInterpolator(const std::array<xType, N>& points, const std::array<yType, N>& values, const std::array<yType, N>& deriv) noexcept {
    for (unsigned int i = 0; i < N; ++i) {
        duplicated_x[2 * i] = points[i];
        duplicated_x[2 * i + 1] = points[i];
        duplicated_f[2 * i] = values[i];
        duplicated_f[2 * i + 1] = values[i];
    }
    for (unsigned int i = 0; i < 2 * N; ++i) {
        divided_diff[i][0] = duplicated_f[i];
    }
    for (unsigned int j = 1; j < 2 * N; ++j) {
        for (unsigned int i = 0; i < 2 * N - j; ++i) {
            if (duplicated_x[i] == duplicated_x[i + j]) {
                divided_diff[i][j] = deriv[i / 2];
            }
            else {
                divided_diff[i][j] = (divided_diff[i + 1][j - 1] - divided_diff[i][j - 1]) / (duplicated_x[i + j] - duplicated_x[i]);
            }
        }
    }
    for (unsigned int j = 0; j < 2 * N; ++j) {
        coefficients[j] = divided_diff[0][j];
    }
}

template<typename xType, typename yType, unsigned int N>
yType NewtonInterpolator<xType, yType, N>::interpolate(const xType& x) const noexcept {
    yType result = coefficients[2 * N - 1];
    for (int j = 2 * N - 2; j >= 0; --j) {
        result = result * (x - duplicated_x[j]) + coefficients[j];
    }
    return result;
}
