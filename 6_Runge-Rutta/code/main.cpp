#include "funcs.h"
#include<fstream>

double func1(double y, double t) { return(t * t * t); }
double ans1(double t) { return (t * t * t * t / 4); }
double func2(double y, double t) { return(-y); }
double ans2(double t) { return (std::sin(t)); }

ÑauchyProblem<decltype(func1), 1> linear;
ÑauchyProblem<decltype(func2), 2> oscillator;

int main() {
    const double start = 0;
    const double end = 5;

    std::ofstream data1("data_I.txt");
    data1.precision(16);
    std::ofstream data2("data_II.txt");
    data2.precision(16);

    for (double h = -7; h < 0; h += 0.1) {
        
        double step = std::pow(10, h);
        
        std::vector<std::array<double, 2>> res1
            = integrate<decltype(func1), RK4Table, ÑauchyProblem<decltype(func1), 1>>(func1, Eigen::Vector<double, 1>(0.0), start, end, step, linear);
        double err1 = std::abs(res1[0][0] - ans1(res1[0][1]));

        double delta;
        for (unsigned int i = 1; i < res1.size(); i++) {

            delta = std::abs(res1[i][0] - ans1(res1[i][1]));
            err1 = delta > err1 ? delta : err1;
        }

        std::vector<std::array<double, 2>> res2
            = integrate<decltype(func2), RK4Table, ÑauchyProblem<decltype(func2), 2>>(func2, Eigen::Vector<double, 2>(0.0, 1.0), start, end, step, oscillator);
        double err2 = std::abs(res2[0][0] - ans2(res2[0][1]));

        for (unsigned int i = 1; i < res2.size(); i++){
            delta = std::abs(res2[i][0] - ans2(res2[i][1]));
            err2 = delta > err2 ? delta : err2;
        }
        
        data1 << h << "     " << err1 << std::endl;
        data2 << h << "     " << err2 << std::endl;
    }

    data1.close();
    data2.close();

    return 0;
};