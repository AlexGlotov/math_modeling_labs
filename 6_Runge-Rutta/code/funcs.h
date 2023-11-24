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

struct DP45 {
    static constexpr unsigned int stages = 7;
    static constexpr std::array<std::array<double, stages>, stages> table = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                                              0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                                              0.075, 0.225, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                                              44.0/45.0, -56.0/15.0, 32.0/9.0, 0.0, 0.0, 0.0, 0.0, 
                                                                              19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0.0, 0.0, 0.0,
                                                                              9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0, 0.0, 0.0,
                                                                              35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0 };
    static constexpr std::array<double, stages> cColumn = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
    static constexpr std::array<double, stages> bString1 = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0};
    static constexpr std::array<double, stages> bString2 = {5179.0/57600, 0.0, 7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 0.025};
    static constexpr unsigned int approximation = 5;
};

struct StepControl {
    double minStep;
    double maxStep;
    double tolerance;
    double initialStep;
};

template<typename Callable, typename Table, typename RHS>  //интегрирование с контролем шага
std::vector<typename RHS::StateAndArg> integrateDP(
    const Callable& func,
    const typename RHS::StateAndArg& initialState,
    const double& endTime,
    const StepControl& stepControl,
    const RHS& rhs)
{

    const unsigned int stages = Table::stages;
    const unsigned int dim = RHS::dim;
    const unsigned int p = Table::approximation;
    const double tol = stepControl.tolerance;

    typename RHS::StateAndArg state = initialState;

    Eigen::Matrix<double, dim, stages> k;

    std::vector<typename RHS::StateAndArg> res;
    res.push_back(state);

    double delta = stepControl.tolerance * 10;
    for (double t = initialState.arg; t < endTime; )
    {
        Eigen::Vector < double, dim> dY2;
        double step = stepControl.initialStep;
        for (step; delta > tol and
            step > stepControl.minStep; )
        {
            step *= std::pow(tol / delta, 1.0 / p);
            //1sf int
            for (int i = 0; i < stages; i++) {
                Eigen::Vector<double, dim> Y = state.state;
                for (int j = 0; j < i; j++) {
                    Y += Table::table[i][j] * k.col(j);
                }
                k.col(i) = rhs.calc(func, { Y, state.arg + step * Table::cColumn[i] });
                k.col(i) = step * k.col(i);
            }

            Eigen::Vector < double, dim> dY1;
            dY1 = Eigen::Vector < double, dim>::Zero();
            for (int i = 0; i < stages; i++) {
                dY1 += Table::bString1[i] * k.col(i);
            }

            //2nd int
            for (int i = 0; i < stages; i++) {
                Eigen::Vector<double, dim> Y = state.state;
                for (int j = 0; j < i; j++) {
                    Y += Table::table[i][j] * k.col(j);
                }
                k.col(i) = rhs.calc(func, { Y, state.arg + step * Table::cColumn[i] });
                k.col(i) = step * k.col(i);
            }

            
            dY2 = Eigen::Vector < double, dim>::Zero();
            for (int i = 0; i < stages; i++) {
                dY2 += Table::bString2[i] * k.col(i);
            }

            delta = (dY1 - dY2).norm();

        }

        t += step;
        state.state += dY2;
        state.arg += step;
        res.push_back(state);
        delta = stepControl.tolerance * 10;

    }

    return res;
};


//y^(n) = f(t, y)
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
    Eigen::Vector<double, dim> calc(const Callable& func, const State& value, const Argument arg) const {
        Eigen::Vector<double, dim> res;
        for (int i = 0; i < dim - 1; i++) {
            res[i] = value[i + 1];
        }
        res[dim - 1] = func(value[0], arg);
        return res;
    };

    Eigen::Vector<double, dim> calc(const Callable& func, const StateAndArg state) const {
        return (func(state.state));

    };
};



template<typename Callable, typename Table, typename RHS>  // интегрирование задачи Коши одной переменной без контроля шага
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
