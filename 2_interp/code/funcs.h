#pragma once

#include <iostream>
#include <array>


template <typename RealType, unsigned int N>
std::array<RealType, N> linear_grid(const double start, const double final) {
    std::array<RealType, N> arr;
    double step = (final - start) / (N - 1);

    arr[0] = start;
    for (unsigned int i = 1; i < N; ++i) {
        arr[i] = arr[i - 1] + step;
    }
    return arr;
};

template <typename RealType, unsigned int N>
std::array<RealType, N> chebyshev_grid(const double start, const double final) {
    std::array<RealType, N> arr;

    const double pi = acos(-1);

    double h_s = (start + final) / 2;
    double h_d = (final - start) / 2;

    for (unsigned int i = 0; i < N; i++) {
        arr[i] = h_s + h_d * cos(pi * (2 * (N - 1 - i) + 1) / (2 * N));
    }
    return arr;
}


template<typename xType, typename yType, unsigned int N>
class HermiteInterpolator {
private:

    std::array<yType, 2 * N> coeffs;
    std::array<xType, N> grid;

public:

    HermiteInterpolator(const std::array<xType, N>& points,
        const std::array<yType, N>& values, const std::array<yType, N>& deriv) noexcept :
        grid{ points }
    {
        coeffs[0] = values[0];

        for (unsigned int i = 1; i < N; i++) {
            coeffs[2 * i] = (values[i] - values[i - 1]) / (grid[i] - grid[i - 1]);
            coeffs[2 * i - 1] = deriv[i - 1];
        }
        coeffs[2 * N - 1] = deriv[N - 1];

        for (unsigned int i = 1; i < 2 * N; i++) {
            for (unsigned int j = 2 * N - 1; j > i; j--) {
                coeffs[j] = (coeffs[j] - coeffs[j - 1]) / (grid[j / 2] - grid[(j - i - 1) / 2]);
            }
        }

    }
    
    yType interpolate(const xType& x) const noexcept
    {
        yType y;

        y = coeffs[2 * N - 1];

        for (int i = 2 * N - 2; i >= 0; i--)
        {
            y = coeffs[i] + y * (x - grid[i / 2]);
        }

        return y;
    };
   
};


template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
private:

    std::array<yType, N> coeffs;
    std::array<xType, N> grid;

public:

    NewtonInterpolator(const std::array<xType, N>& points, const std::array<yType, N>& values) noexcept :
        coeffs{ values }, grid{ points }
    {
        for (unsigned int i = 0; i < N; i++) {
            for (unsigned int j = N - 1; j > i; j--) {
                coeffs[j] = (coeffs[j] - coeffs[j - 1]) / (grid[j] - grid[j - i - 1]);
            }
        }
    };

    yType interpolate(const xType& x) const noexcept
    {
        yType y;

        y = coeffs[N - 1];

        for (int i = N - 2; i >= 0; i--) {
            y = coeffs[i] + y * (x - grid[i]);
        }

        return y;
    };
};
