#include "funcs.h"
#include <fstream>
#include <array>
#include <cmath>

int main() {
    const double referencePoint = 1.0;
    const std::array<double, 16> stepValues = { 
        1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 
        1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 
        1e-11, 1e-12, 1e-13, 1e-14, 1e-15 
    };
    const double exponential = std::exp(1.0);

    std::ofstream file2("file2.txt");
    if (!file2) return -1;
    file2.precision(10);
    for (const auto& step : stepValues) {
        file2 << step << "\t" << d_e<double, 2>(referencePoint, step, { -1, 1 }) << "\n";
    }
    file2.close();

    std::array<double, 100000> logarithmicSteps;
    logarithmicSteps[0] = -35.0;
    const double increment = 0.00035;
    for (size_t i = 1; i < logarithmicSteps.size(); ++i) {
        logarithmicSteps[i] = logarithmicSteps[i - 1] + increment;
    }

    std::ofstream file3("file3.txt");
    if (!file3) return -1;
    file3.precision(10);
    for (const auto& x : logarithmicSteps) {
        double derivative = d_e<double, 3>(referencePoint, std::exp(x), { -1, 1, 2 });
        double error = std::log(std::abs(derivative - exponential));
        file3 << x << " " << error << "\n";
    }
    file3.close();

    std::ofstream file4("file4.txt");
    if (!file4) return -1;
    file4.precision(10);
    for (const auto& x : logarithmicSteps) {
        double derivative = d_e<double, 4>(referencePoint, std::exp(x), { -2, -1, 1, 2 });
        double error = std::log(std::abs(derivative - exponential));
        file4 << x << " " << error << "\n";
    }
    file4.close();

    std::ofstream file5("file5.txt");
    if (!file5) return -1;
    file5.precision(10);
    for (const auto& x : logarithmicSteps) {
        double derivative = d_e<double, 5>(referencePoint, std::exp(x), { -2, -1, 1, 2, 3 });
        double error = std::log(std::abs(derivative - exponential));
        file5 << x << " " << error << "\n";
    }
    file5.close();

    return 0;
}
