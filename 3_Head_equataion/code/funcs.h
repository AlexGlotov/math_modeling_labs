#pragma once
#include<vector>
#include<iostream>
#include<type_traits>
#include<array>

struct cell {
    double x;
    double u;
};

struct mesh1D {
    std::vector<cell> layer;
    double t;
};

template<typename Callable>
std::vector<double> Create_Mesh(double x_start, double x_end, double step, unsigned int N, const Callable u_0) {
    std::vector<double> data;


    double x = 0;
    for (unsigned int i = 0; i <= N; i++) {
        data.push_back(u_0(x));
        x += step;
    }

    return data;
}

struct BoundCond
{
    double a_left;
    double b_left;
    double a_right;
    double b_right;
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

template<typename mType>
class ThreeDiagonalMatrix {
private:
    std::vector<std::array<mType, 3>> diag;
public:
    ThreeDiagonalMatrix(const double step, const BoundCond consts, const double CFL, const int N) {

        double a_left = consts.a_left;
        double b_left = consts.b_left;
        double a_right = consts.a_right;
        double b_right = consts.b_right;

        diag.push_back({ 0, a_left - b_left / step , b_left / step - b_left / (2 * step * CFL)});
        for (unsigned int i = 1; i < N; i++)
        {
            diag.push_back({ - CFL, 1.0 + 2.0 * CFL, - CFL });
        }
        diag.push_back({ - b_right / step + b_right / (2 * step * CFL), a_right + b_right / step, 0});

    }

    std::array<double, 3> operator() (const unsigned int i) const {
        return diag[i];
    };

};

template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType>& matrix,
    const std::vector<cType>& column)
{
    unsigned int N = column.size();
    int n = N - 1;
    std::vector<mType> p_vector{ 0 }, q_vector{ 0 };
    p_vector.reserve(N);
    q_vector.reserve(N);
    for (int i = 0; i < N; ++i) {
        p_vector.push_back(- matrix(i)[2] /
            (matrix(i)[0] * p_vector[i] + matrix(i)[1]));
        q_vector.push_back((column[i] - matrix(i)[0] * q_vector[i]) /
            (matrix(i)[0] * p_vector[i] + matrix(i)[1]));
    }
    std::vector<DivisType<cType, mType>> sol(N);
    sol[n] = (column[n] - matrix(n)[0] * q_vector[n]) /
        (matrix(n)[0] * p_vector[n] + matrix(n)[1]);
    for (int i = n - 1; i >= 0; --i) {
        sol[i] = p_vector[i + 1] * sol[i + 1] + q_vector[i + 1];
    }
    return sol;
};

template<typename Callable_bound, typename Callable_func>
std::vector<double> create_column(std::vector<double> data, const double CFL, double step,
                                   Callable_func func, Callable_bound F_left, Callable_bound F_right,
                                   double x_start, double tau, double t, double b_left, double b_right)
{
    std::vector<double> col;
    double S = data.size();

    col.push_back(0);
    for (int i = 1; i < S - 1; i++)
    {
        col.push_back(CFL * data[i - 1] + CFL * data[i + 1] + (1- 2 * CFL) * data[i] + 
                        tau * (func(x_start + i * step, t + tau) + func(x_start + i * step, t)) / 2);

    }

    col.push_back(F_right(t + tau) + col[S - 2] * b_right / (2 * step * CFL));

    col[0] = F_right(t + tau) - col[1] * b_left / (2 * step * CFL);

    return col;
}