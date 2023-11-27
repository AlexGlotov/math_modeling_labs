#pragma once
#include"Eigen/Dense"
#include<vector>
#include<iostream>
#include<array>

struct RK4Table {
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = { 0, 0, 0, 0,
                                                                              0.5, 0, 0, 0,
                                                                              0, 0.5, 0, 0,
                                                                              0, 0, 1.0, 0 };
    static constexpr std::array<double, stages> cColumn = { 0, 0.5, 0.5, 1 };
    static constexpr std::array<double, stages> bString = { double(1) / 6, double(1) / 3, double(1) / 3, double(1) / 6 };
};

struct BDF4 {
    static constexpr unsigned int size = 4;
    static constexpr std::array<double, size + 1> alpha = {1.92, -1.44, 0.64, -0.12, 0.48};
    static constexpr std::array<double, size> beta = { 55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -9.0 / 24.0 };
};

template<typename Callable, unsigned int d>
class СauchyProblem {
public:

    static constexpr unsigned int dim = d;

    using State = Eigen::Vector<double, dim>;
    using Argument = double;

    struct StateAndArg {
        State state;
        Argument arg;
    };

    /* Вычисляет правую часть ДУ - функцию f*/
    Eigen::Vector<double, dim> calc(const Callable& func, const StateAndArg& state) const {
        return (func(state.state, state.arg));
    };

    /*** Вычисляет разницу двух состояний (решения нелиненого уравнения) ***/
    double calcDif(const State& first, const State& second) const {
        return ((second - first).norm());
    };
};

struct IntegrationParameters {
    double step;  // шаг интегрирования
    double epsilon;  // точность решения нелинейного уравнения
    unsigned int maxIter;  // максимальное количество итераций для решения нелинейного уравнения
};

/***
BDF - структура с коэффициентами метода
RHS - правая часть Д.У.
RKTable -  таблица Бутчера метода Рунге-Кутты для разгона.
***/

template<typename Callable, typename BDF, typename RHS, typename RKTable>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
    const Callable& func,
    const typename RHS::StateAndArg& initialState,
    const typename RHS::Argument& endTime,
    const IntegrationParameters& parameters,
    const RHS& rhs)
{
    const unsigned int dim = RHS::dim;
    const unsigned int stages = RKTable::stages;
    const unsigned int size = BDF::size;

    typename RHS::StateAndArg state = initialState;

    const unsigned int N = (endTime - initialState.arg) / parameters.step + 1;
    const double h = (endTime - initialState.arg) / N;
    const unsigned int maxIter = parameters.maxIter;
    const double tol = parameters.epsilon;

    Eigen::Matrix<double, dim, stages> k;

    std::vector<typename RHS::StateAndArg> res;
    res.push_back(state);

    for (unsigned int m = 0; m < size + 1; m++) {

        for (int i = 0; i < stages; i++) {
            Eigen::Vector<double, dim> Y = state.state;
            for (int j = 0; j < i; j++) {
                Y += RKTable::table[i][j] * k.col(j);
            }
            k.col(i) = rhs.calc(func, { Y, state.arg + h * RKTable::cColumn[i] });
            k.col(i) = h * k.col(i);
        }

        Eigen::Vector < double, dim> dY;
        dY = Eigen::Vector < double, dim>::Zero();
        for (int i = 0; i < stages; i++) {
            dY += RKTable::bString[i] * k.col(i);
        }

        state.state += dY;
        state.arg += h;
        res.push_back(state);
    }
    
    const std::array<double, size + 1> alpha = BDF::alpha;
    const std::array<double, size> beta = BDF::beta;

    for (unsigned int m = size + 1; m < N; m++) {
        Eigen::Vector<double, dim> Y1 = state.state;

        for (int l = 0; l < size; l++) {
            Y1 += h * beta[l] * func(res[m - l].state, res[m - l].arg);
        }

        double delta = tol * 10;
        Eigen::Vector<double, dim> Y2;

        for (int i = 0; i < maxIter and delta > tol; i++) {
            Y2 = h * alpha[size] * func(Y1, state.arg + h);
            for (int j = 0; j < size; j++) {
                Y2 += alpha[j] * res[m - j].state;
            }

            delta = rhs.calcDif(Y1, Y2);

            Y1 = Y2;
        }

        state.state = Y2;
        state.arg += h;

        res.push_back(state);

    }

    return res;
};