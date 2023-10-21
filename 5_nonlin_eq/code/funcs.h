#pragma once

#include<iostream>
#include <array>
#include <type_traits>
#include <tuple>

struct IterationException {};

/**
    Решает уравнение Кеплера методом Ньютона
    * ecc - эксцентриситет, принадлежит (0, 1)
    * meanAnomaly - средняя аномалия, М (радианы)
    * maxIter - максимальное количество итераций
    * tol - точность, с которой нужно отыскать решение

    Рекомендуемое поведение. Если решение не нашлось за maxIter итераций - выбрасывать исключение.
    Если приближения к решению между итерациями меняются не более, чем на tol, то решение достигнуто.
**/
double keplerSolver(double ecc, double meanAnomaly, unsigned int maxIter, double tol)
{
    double E = meanAnomaly + ecc * std::sin(meanAnomaly) /
        (1 - ecc * std::cos(meanAnomaly));
    double delta = std::abs(meanAnomaly - E);

    unsigned int i;
    for (i = 1; i < maxIter and delta > tol; i++) {
        delta = E;
        E = E - (E - ecc * std::sin(E) - meanAnomaly) /
            (1 - ecc * std::cos(E));
        delta = std::abs(delta - E);
    }
    
    //if (i == maxIter) { throw IterationException(); };

    return E;
};



template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename Callable, typename RealType>
decltype(auto) solve(
    const Callable& func,                                             // функция F
    const RealType& tau,                                              // шаг тау
    const typename ArgumentGetter<Callable>::Argument& initialGuess,  // начальное приближение
    const unsigned int nIteration                                     // количество итераций
) 
{
    double x = initialGuess;

    for (int i = 0; i < nIteration; i++) {
        x += tau * func(x);
    }

    return x;

};