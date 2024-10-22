#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cassert>

// Structure to hold derivative coefficients
template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoeffs;
};

// Function to solve linear systems using Gaussian elimination
template<typename RealType>
bool solveLinearSystem(std::vector<std::vector<RealType>>& A, std::vector<RealType>& b, std::vector<RealType>& solution) {
    const unsigned int n = A.size();
    for(unsigned int i = 0; i < n; ++i){
        // Partial pivoting
        RealType maxElem = std::abs(A[i][i]);
        unsigned int maxRow = i;
        for(unsigned int k = i+1; k < n; ++k){
            if(std::abs(A[k][i]) > maxElem){
                maxElem = std::abs(A[k][i]);
                maxRow = k;
            }
        }
        if(maxElem < 1e-12){
            // Singular matrix
            return false;
        }
        // Swap maximum row with current row (pivot)
        std::swap(A[i], A[maxRow]);
        std::swap(b[i], b[maxRow]);

        // Eliminate below
        for(unsigned int k = i+1; k < n; ++k){
            RealType factor = A[k][i] / A[i][i];
            for(unsigned int j = i; j < n; ++j){
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Back substitution
    solution.assign(n, 0);
    for(int i = n-1; i >=0; --i){
        RealType sum = 0;
        for(unsigned int j = i+1; j < n; ++j){
            sum += A[i][j] * solution[j];
        }
        solution[i] = (b[i] - sum) / A[i][i];
    }
    return true;
}

// Function to calculate derivative coefficients using the method of undetermined coefficients
template<typename RealType, unsigned int N, unsigned int L>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {
    constexpr unsigned int equations = N + 1;
    std::vector<std::vector<RealType>> A(equations, std::vector<RealType>(N + 1, 0));
    std::vector<RealType> b(equations, 0);

    // Set up the system
    for(unsigned int m = 0; m < equations; ++m){
        // Central node (k=0)
        A[m][0] = std::pow(0, m); // which is 1 if m=0, else 0
        // Other nodes
        for(unsigned int i = 0; i < N; ++i){
            A[m][i+1] = std::pow(points[i], m);
        }
        // Right-hand side
        if(m == L){
            b[m] = std::tgamma(L + 1); // L! using gamma function
        } else {
            b[m] = 0;
        }
    }

    // Solve the linear system
    std::vector<RealType> solution;
    bool success = solveLinearSystem(A, b, solution);
    assert(success && "Linear system could not be solved.");

    DerivativeCoef<RealType, N> coef;
    coef.centralCoef = solution[0];
    for(unsigned int i = 0; i < N; ++i){
        coef.otherCoeffs[i] = solution[i+1];
    }
    return coef;
}

// Function to generate points based on N (symmetric if possible)
std::vector<double> generatePoints(unsigned int N){
    std::vector<double> points;
    if(N == 0){
        return points;
    }
    unsigned int half = N / 2;
    for(int i = half; i >=1; --i){
        points.push_back(-static_cast<double>(i));
    }
    if(N %2 !=0){
        points.push_back(-0.5);
    }
    for(unsigned int i =1; i <= half; ++i){
        points.push_back(static_cast<double>(i));
    }
    if(N %2 !=0){
        points.push_back(0.5);
    }
    return points;
}

int main(){
    // Define the exact second derivative of y = e^x at x=1
    const double exactDerivative = std::exp(1.0);

    // Define step sizes h from 1 to 1e-15
    std::vector<double> hs;
    for(int i =0; i <=15; ++i){
        hs.push_back(std::pow(10.0, -static_cast<double>(i)));
    }

    // Define N values
    std::vector<unsigned int> Ns = {3,4,5};

    // Set precision for output
    std::cout << std::fixed << std::setprecision(8);

    // Header
    std::cout << "N, h, Approximation, Error, log(h), log(Error)\n";

    // Storage for slopes
    std::vector<std::pair<unsigned int, std::vector<std::pair<double, double>>>> data;

    for(auto N : Ns){
        // Generate points
        std::vector<double> points = generatePoints(N);
        std::array<double, 5> arrPoints;
        for(unsigned int i =0; i < N && i < arrPoints.size(); ++i){
            arrPoints[i] = points[i];
        }

        // Calculate derivative coefficients for L=2
        DerivativeCoef<double, 5> coef;
        if(N ==2){
            std::array<double,2> pts = {points[0], points[1]};
            coef = calcDerivativeCoef<double,2,2>(pts);
        }
        else if(N ==3){
            std::array<double,3> pts = {points[0], points[1], points[2]};
            coef = calcDerivativeCoef<double,3,2>(pts);
        }
        else if(N ==4){
            std::array<double,4> pts = {points[0], points[1], points[2], points[3]};
            coef = calcDerivativeCoef<double,4,2>(pts);
        }
        else if(N ==5){
            std::array<double,5> pts = {points[0], points[1], points[2], points[3], points[4]};
            coef = calcDerivativeCoef<double,5,2>(pts);
        }

        // Store log(h) and log(error)
        std::vector<std::pair<double, double>> logData;

        for(auto h : hs){
            // Define function values at the nodes
            // f(x0 + k_i h) = e^(1 + k_i h)
            std::vector<double> fValues;
            for(unsigned int i =0; i < N; ++i){
                fValues.push_back(std::exp(1.0 + points[i]*h));
            }
            double f0 = std::exp(1.0); // f(x0)

            // Approximate derivative: (centralCoef * f0 + sum(otherCoeffs * f(x0 +k_i h))) / h^2
            double approximation = (coef.centralCoef * f0);
            for(unsigned int i =0; i < N; ++i){
                approximation += coef.otherCoeffs[i] * fValues[i];
            }
            approximation /= (h*h);

            // Compute error
            double error = std::abs(approximation - exactDerivative);

            // Compute log(h) and log(error)
            double logh = std::log10(h);
            double logError = (error > 0) ? std::log10(error) : -100; // Assign a large negative value if error is zero

            // Output the data
            std::cout << N << ", " << h << ", " << approximation << ", " << error << ", " << logh << ", " << logError << "\n";

            // Store for slope calculation
            logData.emplace_back(logh, logError);
        }
        data.emplace_back(N, logData);
    }

    // Note: Plotting and slope estimation would typically be done using external tools like Python's matplotlib.
    // Below is a simple slope estimation using linear regression on the log-log data for each N.

    std::cout << "\nEstimated Slopes for Different N:\n";
    for(auto& entry : data){
        unsigned int N = entry.first;
        auto& logData = entry.second;

        // To estimate slope, perform linear regression on logh vs logError
        // Considering only the linear region. For simplicity, we'll consider the last 5 points (small h)
        size_t m = logData.size();
        size_t start = m >5 ? m -5 : 0;
        double sumx=0, sumy=0, sumxy=0, sumx2=0;
        size_t count =0;
        for(size_t i = start; i < m; ++i){
            double x = logData[i].first;
            double y = logData[i].second;
            if(y < -100) continue; // Ignore invalid log errors
            sumx += x;
            sumy += y;
            sumxy += x * y;
            sumx2 += x * x;
            count++;
        }
        if(count <2){
            std::cout << "N=" << N << ": Not enough data for slope estimation.\n";
            continue;
        }
        double slope = (count * sumxy - sumx * sumy) / (count * sumx2 - sumx * sumx);
        std::cout << "N=" << N << ": Slope ≈ " << slope << "\n";
    }

    return 0;
}
