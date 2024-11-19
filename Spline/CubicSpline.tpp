#ifndef CUBICSPLINE_TPP
#define CUBICSPLINE_TPP

#include "CubicSpline.hpp"

/** Реализация функции прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType>& matrix,
                                           const std::vector<cType>& column) {
    int n = matrix.main.size();
    std::vector<DivisType<cType, mType>> c_prime(n, 0);
    std::vector<DivisType<cType, mType>> d_prime(n, 0);

    if (matrix.main[0] == 0) {
        throw std::runtime_error("Zero on main diagonal");
    }

    c_prime[0] = matrix.upper[0] / matrix.main[0];
    d_prime[0] = column[0] / matrix.main[0];

    for(int i = 1; i < n; ++i){
        mType denom = matrix.main[i] - matrix.lower[i-1] * c_prime[i-1];
        if (denom == 0) {
            throw std::runtime_error("Zero pivot encountered");
        }
        if(i < n-1){
            c_prime[i] = matrix.upper[i] / denom;
        }
        d_prime[i] = (column[i] - matrix.lower[i-1] * d_prime[i-1]) / denom;
    }

    std::vector<DivisType<cType, mType>> x(n, 0);
    x[n-1] = d_prime[n-1];
    for(int i = n-2; i >=0; --i){
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }

    return x;
}

/** Реализация конструктора CubicSpline **/
template<typename xType, typename yType>
CubicSpline<xType, yType>::CubicSpline(const std::vector<xType>& pts, 
                                      const std::vector<yType>& vals, 
                                      const Deriv2Type& first, 
                                      const Deriv2Type& second)
    : points(pts), values(vals)
{
    if(points.size() != values.size()){
        throw std::invalid_argument("Points and values must have the same size.");
    }

    int n = points.size();
    if(n < 2){
        throw std::invalid_argument("At least two points are required.");
    }

    std::vector<DeltaXType> h(n-1);
    for(int i = 0; i < n-1; ++i){
        h[i] = points[i+1] - points[i];
        if(h[i] <= 0){
            throw std::invalid_argument("Points must be in strictly increasing order.");
        }
    }

    std::vector<DerivType> alpha(n, 0.0);
    for(int i =1; i < n-1; ++i){
        alpha[i] = (3.0 * (values[i+1] - values[i])/h[i] - 
                    3.0 * (values[i] - values[i-1])/h[i-1]);
    }

    ThreeDiagonalMatrix<Deriv2Type> A(n);
    // Заполнение диагоналей
    A.main[0] = 1.0;
    A.upper[0] = 0.0;
    A.lower[n-2] = 0.0;
    A.main[n-1] = 1.0;

    for(int i =1; i < n-1; ++i){
        A.lower[i-1] = h[i-1];
        A.main[i] = 2.0*(h[i-1] + h[i]);
        A.upper[i] = h[i];
    }

    std::vector<DerivType> b(n, 0.0);
    for(int i =1; i < n-1; ++i){
        b[i] = alpha[i];
    }
    // Граничные условия
    b[0] = first;
    b[n-1] = second;

    secondDerivatives = solve(A, b);
}

/** Реализация метода interpolate **/
template<typename xType, typename yType>
yType CubicSpline<xType, yType>::interpolate(const xType& x) const noexcept {
    int n = points.size();
    if(x < points.front() || x > points.back()){
        return NAN; // Или можно выбросить исключение
    }

    // Поиск интервала с помощью бинарного поиска
    int low = 0;
    int high = n -1;
    int mid = 0;
    while(low <= high){
        mid = low + (high - low) / 2;
        if(mid == n-1){
            break;
        }
        if(points[mid] <= x && x < points[mid+1]){
            break;
        }
        if(x < points[mid]){
            high = mid -1;
        }
        else{
            low = mid +1;
        }
    }

    if(mid == n-1){
        mid = n-2;
    }

    DeltaXType h = points[mid+1] - points[mid];
    DerivType a = (points[mid+1] - x) / h;
    DerivType b = (x - points[mid]) / h;
    yType y = a * values[mid] + b * values[mid+1] +
             ((a*a*a - a) * secondDerivatives[mid] + 
              (b*b*b - b) * secondDerivatives[mid+1]) * (h*h) / 6.0;
    return y;
}

#endif // CUBICSPLINE_TPP
