#pragma once
#include <vector>
#include<iostream>
#include <type_traits>
#include<array>

template <typename RealType>
std::vector<RealType> linear_grid(const RealType start, const RealType final, const unsigned int N) {
    std::vector<RealType> arr;
    RealType step = (final - start) / (N - 1);

    arr.push_back(start);
    for (unsigned int i = 1; i < N; ++i) {
        arr.push_back(arr[i - 1] + step);
    }
    return arr;
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

/** класс для работы с трехдиагональной матрицей **/
template<typename Type>
class ThreeDiagonalMatrix {
private:
    std::vector<std::array<Type, 3>> diag;
public:
    ThreeDiagonalMatrix(const std::vector<Type> &step) {
        const unsigned int N = step.size();
        
        diag.reserve(N - 1);

        diag.push_back({ 0, 2.0 / 3 , step[1] / (3 * (step[0] + step[1])) });
        for (unsigned int i = 1; i < N - 2; i++)
        {
            diag.push_back({ step[i + 1] / (3 * (step[i] + step[i + 1])),
                2.0 / 3, step[i] / (3 * (step[i] + step[i + 1])) });
        }
        diag.push_back({ step[N - 2] / (3 * (step[N - 2] + step[N - 1])), 2.0 / 3, 0 });
    
    }

    //вызов i строки 
    std::array<Type, 3> operator() (const unsigned int i) const {
        return diag[i];
    };

};

/** Функция для решения метод0м прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType>& matrix,
    const std::vector<cType>& column) 
{
    const unsigned int N = column.size();

    std::vector<mType> p, q;
    p.reserve(N - 1); q.reserve(N - 1);

    p.push_back(-matrix(0)[2] / matrix(0)[1]);
    q.push_back(column[0] / matrix(0)[1]);

    std::vector<DivisType<cType, mType>> sol(N + 2);

    for (unsigned int i = 1; i < N - 1; i++)
    {
        p.push_back(-matrix(i)[2] / (matrix(i)[0] * p[i - 1] + matrix(i)[1]));
        q.push_back((column[i] - matrix(i)[0] * q[i - 1]) / (matrix(i)[0] * p[i - 1] + matrix(i)[1]));
    }

    sol[0] = 0;
    sol[N + 1] = 0;
    sol[N] = (column[N - 1] - matrix(N - 1)[0] * q[N - 2]) /
        (matrix(N - 1)[0] * p[N - 2] + matrix(N - 1)[1]);
    for (unsigned int i = N - 1; i > 0; i--) {
        sol[i] = sol[i + 1] * p[i - 1] + q[i - 1];
    }

    return sol;
};

template<typename xType, typename yType>
std::vector<yType> devided_diffrences(const std::vector<xType>& h,
    const std::vector<yType>& values) {
   
    const unsigned int N = h.size();

    std::vector<yType> dev_dif;
    dev_dif.reserve(N - 1);

    for (unsigned int i = 1; i < N; i++)
    {
        dev_dif.push_back((values[i + 1] * h[i - 1] - values[i] * (h[i] + h[i - 1]) + values[i - 1] * h[i])
            / (h[i] * h[i - 1] * (h[i] + h[i - 1])));
    }

    return dev_dif;
}

template<typename xType, typename yType>
class CubicSpline {
    using DeltaXType = DiffType<xType>;
    using DerivType = DivisType<DiffType<yType>, DeltaXType>;
    using Deriv2Type = DivisType<DiffType<DerivType>, DeltaXType>;
private:
    std::vector<std::array<yType, 4>> coefs;
    std::vector<xType> points;

public:
    CubicSpline(const std::vector<xType>& points, const std::vector<yType>& values,
        const Deriv2Type& first,  // значение для левой второй производной
        const Deriv2Type& second  // значение для правой второй производной
    ) : points{ points }
    {
        const unsigned int N = points.size() - 1;
        coefs.reserve(N);

        std::vector<xType> h;
        h.reserve(N);

        for (unsigned int i = 0; i < N; i++) {
            h.push_back(points[i + 1] - points[i]);
        }

        std::vector<yType> column = devided_diffrences<double, double>(h, values);
        column[0] = column[0] - first * h[0] / (6 * (h[0] + h[1]));
        column[N - 2] = column[N - 2] - second * h[N - 1] / (6 * (h[N - 1] + h[N - 2]));
        std::vector<yType> vec_c = solve(ThreeDiagonalMatrix<double>(h), column);
        vec_c[0] = first / 2;
        vec_c[N] = second / 2;

        for (unsigned int i = 0; i < N; i++) {
            coefs.push_back({ values[i + 1],
                2 * vec_c[i + 1] * h[i] / 3 + vec_c[i] * h[i] / 3 + (values[i + 1] - values[i]) / h[i],
                vec_c[i + 1],
                (vec_c[i + 1] - vec_c[i]) / (3 * h[i]) });
        }

    };

    yType interpolate(const xType& x) const noexcept {

        unsigned int k = std::upper_bound(points.begin(), points.end(), x) - points.begin();

        return(coefs[k - 1][0] + coefs[k - 1][1] * (x - points[k]) + coefs[k - 1][2] * (x - points[k]) * (x - points[k]) +
            coefs[k - 1][3] * (x - points[k]) * (x - points[k]) * (x - points[k]));
    };
};