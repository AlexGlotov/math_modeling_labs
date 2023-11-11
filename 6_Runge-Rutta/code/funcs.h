#pragma once

#include<Eigen\Dense>
#include<vector>
#include<iostream>
#include<array>



/*таблица Бутчера для метода Рунге-Кутты 4 порядка*/
struct RK4Table {
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = { 0, 0, 0, 0, 
                                                                              0.5, 0, 0, 0, 
                                                                              0, 0.5, 0, 0, 
                                                                              0, 0, 1.0, 0 };
    static constexpr std::array<double, stages> cColumn = { 0, 0.5, 0.5, 1 };
    static constexpr std::array<double, stages> bString = {double(1)/6, double(1)/3, double(1)/3, double(1)/6 };
};

//y^(n) = f(t, y)
template<typename Callable, unsigned int d>
class СauchyProblem {
public:

    static constexpr unsigned int dim = d;

    using State = Eigen::Vector<double, dim>;
    using Argument = double;

    State state; //variables
    Argument arg; //argument

    /* Вычисляет правую часть ДУ - функцию f*/
    Eigen::Vector<double, dim> calc(const Callable& func, const State& value, const Argument arg) const {
        Eigen::Vector<double, dim> res;
        for (int i = 0; i < dim - 1; i++) {
            res[i] = value[i + 1];
        }
        res[dim - 1] = func(value[0], arg);
        return res;

    };
};

template<typename Callable, typename Table, typename RHS>  // таблица бутчера и класс правой части f
std::vector<std::array<double, 2>> integrate(
    const Callable& func,
    const typename RHS::State& initialState,
    const typename RHS::Argument& initialArg,
    const typename RHS::Argument& endTime,
    double step,
    const RHS& rhs)
{
    const unsigned int dim = RHS::dim;
    const unsigned int stages = Table::stages;

    Eigen::Vector<double, dim> state = initialState;
    double arg = initialArg;
    Eigen::Matrix<double, dim, stages> k;

    const unsigned int N = (endTime - initialArg) / step + 1;
    step = (endTime - initialArg) / N;

    std::vector<std::array<double, 2>> res;
    res.push_back({ state[0], arg });

    for (unsigned int n = 0; n < N; n++) {

        for (int i = 0; i < stages; i++) {
            Eigen::Vector<double, dim> Y = state;
            for (int j = 0; j < i; j++) {
                Y += Table::table[i][j] * k.col(j);
            }
            k.col(i) = rhs.calc(func, Y, arg + step * Table::cColumn[i]);
            k.col(i) = step * k.col(i);
        }

        Eigen::Vector < double, dim> dY;
        dY = Eigen::Vector < double, dim>::Zero();
        for (int i = 0; i < stages; i++) {
            dY += Table::bString[i] * k.col(i);
        }

        state += dY;
        arg += step;
        res.push_back({ state[0], arg });
    }

    return res;
}