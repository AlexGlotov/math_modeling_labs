#pragma once
#include <vector>
#include<iostream>
#include <type_traits>

template <typename RealType>
std::vector<RealType> linear_grid(const RealType start, const RealType final, const unsigned int N) {
    std::vector<RealType> arr;
    double step = (final - start) / (N - 1);

    arr.push_back(start);
    for (unsigned int i = 1; i < N; ++i) {
        arr.push_back(arr[i - 1] + step);
    }
    return arr;
};

/** класс для работы с трехдиагональной матрицей **/
template<typename Type>
class ThreeDiagonalMatrix {
private:
    std::vector<Type> up_vec; 
    std::vector<Type> centr_vec; 
    std::vector<Type> und_vec; 
public:
    ThreeDiagonalMatrix(const std::vector<Type>& a, const std::vector<Type>& b, 
        const std::vector<Type>& c) : up_vec{ a }, centr_vec{ b }, und_vec{ c } {};

    ThreeDiagonalMatrix(const std::vector<Type> step) {
        const unsigned int N = step.size();
        
        centr_vec.push_back(2.0 / 3);
        up_vec.push_back(step[1] / (3 * (step[0] + step[1])));
        for (unsigned int i = 1; i < N - 2; i++)
        {
            centr_vec.push_back(2.0 / 3);
            up_vec.push_back(step[i + 1] / (3 * (step[i] + step[i + 1])));
            und_vec.push_back(step[i] / (3 * (step[i] + step[i + 1])));
        }
        centr_vec.push_back(2.0 / 3);
        und_vec.push_back(step[N - 2] / (3 * (step[N - 2] + step[N - 1])));
    
    }

    Type operator() (const unsigned int i, const unsigned int j) const {
        if (i - j == 0) {
            return centr_vec[i];
        }
        else if (i - j == 1) {
            return up_vec[j];
        }
        else if (j - i == 1) {
            return und_vec[i];
        }
        else {
            return -1; //Сделать try catch на ошибки
        };
    };
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

/** Функция для решения метод0м прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType>& matrix,
    const std::vector<cType>& column) 
{
    const unsigned int N = column.size();

    std::vector<mType> p{ -matrix(0, 1) / matrix(0,0) }, q{ column[0] / matrix(0, 0) };

    std::vector<DivisType<cType, mType>> sol(N + 2);

    for (unsigned int i = 1; i < N - 1; i++)
    {
        p.push_back(-matrix(i, i + 1) / (matrix(i, i - 1) * p[i - 1] + matrix(i, i)));
        q.push_back((column[i] - matrix(i, i - 1) * q[i - 1]) / (matrix(i, i - 1) * p[i - 1] + matrix(i, i)));
    }

    sol[0] = 0;
    sol[N + 1] = 0;
    sol[N] = (column[N - 1] - matrix(N - 1, N - 2) * q[N - 2]) /
        (matrix(N - 1, N - 2) * p[N - 2] + matrix(N - 1, N - 1));
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

    for (unsigned int i = 1; i < N; i++)
    {
        dev_dif.push_back((values[i + 1] * h[i - 1] - values[i] * (h[i] + h[i - 1]) + values[i - 1] * h[i])
            / (h[i] * h[i - 1] * (h[i] + h[i - 1])));
    }

    return dev_dif;
}

template<typename xType, typename yType>
class CubicSpline {
private:
    std::vector<yType> vec_a;
    std::vector<yType> vec_b;
    std::vector<yType> vec_c;
    std::vector<yType> vec_d;
    std::vector<xType> points;

public:
    CubicSpline(const std::vector<xType>& points, const std::vector<yType>& values) : points{ points }
    {
        const unsigned int N = points.size() - 1;

        std::vector<xType>h;

        for (unsigned int i = 0; i < points.size() - 1; i++) {
            h.push_back(points[i + 1] - points[i]);
        }

        vec_c = solve(ThreeDiagonalMatrix<double>(h),
            devided_diffrences<double, double>(h, values));
        
        for (unsigned int i = 0; i < N; i++) {
            vec_d.push_back((vec_c[i + 1] - vec_c[i]) / (3 * h[i]));
            vec_a.push_back(values[i + 1]);
            vec_b.push_back(2 * vec_c[i + 1] * h[i] / 3 + vec_c[i] * h[i] / 3 +
                (values[i + 1] - values[i]) / h[i]);
        }
        vec_c.erase(vec_c.begin());

    };

    yType interpolate(const xType& x) const noexcept {
    
        yType y;

        if ((x < (*(points.begin() - 0.0001))) || (x > (*(points.end() - 1) + 0.0001))) {
            std::cout << "The value is outside the area under consideration" << std::endl;
            system("pause");
            return -1;
        }
            
        unsigned int k = 1;
        while (x > (points[k] + 0.0001))
        {
            k++;
        }

        y = vec_a[k - 1] + vec_b[k - 1] * (x - points[k]) + vec_c[k - 1] * (x - points[k]) * (x - points[k]) +
            vec_d[k - 1] * (x - points[k]) * (x - points[k]) * (x - points[k]);

        return y;

    };
};