#pragma once

#include<iostream>
#include <array>
#include <type_traits>
#include <tuple>

/*
template<typename RealType, unsigned int N>
const std::array<std::array<RealType, 2>, N> points2weights
    (const std::array<RealType, N> points)
{
    std::array<std::array<RealType, 2>, N> point_weight;
    std::array<RealType, N> coeffs;

    for (int k = 0; k < N; k++) {

        coeffs[0] = 1.0;
        for (int i = 1; i < N; i++) { coeffs[i] = 0; }
        double A = 1;

        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j > 0; j--) {
                coeffs[j] = coeffs[j - 1] - coeffs[j] * points[i];
            }
            coeffs[0] = -coeffs[0] * points[i];
            A = A * (points[k] - points[i]);
        }
        for (int i = k + 1; i < N; i++) {
            for (int j = i; j > 0; j--) {
                coeffs[j] = coeffs[j - 1] - coeffs[j] * points[i];
            }
            coeffs[0] = -coeffs[0] * points[i];
            A = A * (points[k] - points[i]);
        }

        double w = 0;
        for (int i = 0; i < N; i = i + 2) {
            w += 2 * coeffs[i] / (i + 1);
        }
        
        point_weight[k] = { points[k], w / A };
    }
    return point_weight;
}
*/

template <typename RealType, unsigned int N> struct nodes;

template <typename RealType> struct nodes<RealType, 1> {
    static constexpr std::array<RealType, 1> p{ 0.0 };
    static constexpr std::array<RealType, 1> w{ 2.0 };
};

template <typename RealType> struct nodes<RealType, 2> {
    static constexpr std::array<RealType, 2> p{ -0.5773502692, 0.5773502692 };
    static constexpr std::array<RealType, 2> w{ 1.0, 1.0 };
};

template <typename RealType> struct nodes<RealType, 3> {
    static constexpr std::array<RealType, 3> p{ -0.7745966692, 0, 0.7745966692 };
    static constexpr std::array<RealType, 3> w{ 0.5555555556, 0.8888888889, 0.5555555556 };
};

template <typename RealType> struct nodes<RealType, 4> {
    static constexpr std::array<RealType, 4> p{ -0.8611363116, -0.3399810436,
                                               0.3399810436, 0.8611363116 };
    static constexpr std::array<RealType, 4> w{ 0.3478548451,0.6521451549,
                                                0.6521451549, 0.3478548451 };
};

template <typename RealType> struct nodes<RealType, 5> {
    static constexpr std::array<RealType, 5> p{ -0.9061798459, -0.5384693101, 0.0,
                                               0.5384693101, 0.9061798459 };
    static constexpr std::array<RealType, 5> w{ 0.2369268851, 0.4786286705, 0.5688888889,
                                                0.4786286705, 0.2369268851 };
};


template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

/* Ôóíêöèÿ ïðîèçâîäèò èíòåãðèðîâàíèå íà îäíîì îòðåçêå */
template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate(
    const Callable& func,  // Èíòåãðèðóåìàÿ ôóíêöèÿ
    const typename ArgumentGetter<Callable>::Argument& start,  // íà÷àëî îòðåçêà
    const typename ArgumentGetter<Callable>::Argument& end,  // êîíåö îòðåçêà
    const std::array<RealType, N>& points,  // Óçëû êâàäðàòóðû íà îòðåçêå [-1, 1]
    const std::array<RealType, N>& weights // Âåñà óçëîâ êâàäðàòóðû íà îòðåçêà [-1, 1]
)
{
    RealType Int = 0;

    RealType semi_dif = (end - start) / 2;
    RealType semi_sum = (end + start) / 2;

    for (int i = 0; i < N; i++) {
        Int += weights[i] * func(semi_sum + semi_dif * points[i]);
    }

    return Int * semi_dif;
};

/* Ôóíêöèÿ ïðîèçâîäèò èíòåãðèðîâàíèå, ðàçáèâàÿ îòðåçîê íà ïîäîòðåçêè äëèíîé íå áîëåå dx */
template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate(
    const Callable& func,  // Èíòåãðèðóåìàÿ ôóíêöèÿ
    const typename ArgumentGetter<Callable>::Argument& start,  // íà÷àëî îòðåçêà
    const typename ArgumentGetter<Callable>::Argument& end,  // êîíåö îòðåçêà
    const std::array<RealType, N>& points,  // Óçëû êâàäðàòóðû íà îòðåçêå [-1, 1]
    const std::array<RealType, N>& weights, // Âåñà óçëîâ êâàäðàòóðû íà îòðåçêà [-1, 1]
    const Dif<typename ArgumentGetter<Callable>::Argument>& dx  // Äëèíà ïîäîòðåçêà
)
{
    unsigned int L = (end - start) / dx + 1;
    double step = (end - start) / L;

    RealType Int = 0;

    for (int i = 0; i < L; i++) {
        Int += integrate<Callable, RealType, N>(func, start + i * step,
            start + (i + 1) * step, points, weights);
    }

    return Int;
};

template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrateRR(
    const Callable& func,  // Èíòåãðèðóåìàÿ ôóíêöèÿ
    const typename ArgumentGetter<Callable>::Argument& start,  // íà÷àëî îòðåçêà
    const typename ArgumentGetter<Callable>::Argument& end,  // êîíåö îòðåçêà
    const std::array<RealType, N>& points,  // Óçëû êâàäðàòóðû íà îòðåçêå [-1, 1]
    const std::array<RealType, N>& weights, // Âåñà óçëîâ êâàäðàòóðû íà îòðåçêà [-1, 1]
    const RealType err  // âåëè÷èíà îøèáêè
)
{
    double h = end - start;

    RealType Ih = integrate<Callable, RealType, N>(func, start, end, points, weights, h);
    RealType Ih_2 = integrate<Callable, RealType, N>(func, start, end, points, weights, h/2);
    RealType delta = h;
    RealType delta_2 = std::abs(Ih - Ih_2);

    for (double step = h / 4; delta > err; step = step / 2) {
        delta = delta_2;
        Ih = Ih_2;
        Ih_2 = integrate<Callable, RealType, N>(func, start, end, points, weights, step);
        delta_2 = std::abs(Ih - Ih_2);
   }

    double p = std::log(delta / delta_2) / std::log(2);

    return Ih_2 + (Ih - Ih_2) / (pow(2, p) - 1) ;
};

//bonus task
template<typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate_(
    const Callable& func,  // Интегрируемая функция
    const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
    const typename ArgumentGetter<Callable>::Argument& end,  // конец отрезка
    const std::array<RealType, N>& points,  // Узлы квадратуры на отрезке [-1, 1]
    const std::array<RealType, N>& weights, // Веса узлов квадратуры на отрезка [-1, 1]
    const Dif<typename ArgumentGetter<Callable>::Argument>& dx  // Длина подотрезка
)
{
    unsigned int L = (end - start) / dx + 1;
    double step = (end - start) / L;

    unsigned int num_sections = 0;

    RealType Int = 0;

    for (int i = 0; i < L; i++) {
        Int += integrate<Callable, RealType, N>(func, start + i * step,
            start + (i + 1) * step, points, weights);
        num_sections++;
    }
    
    std::tuple<RealType, unsigned int> ret(Int, num_sections);

    return ret;
};
