#include "funcs.h"
#include<fstream>

double func1(double y, double t) { return(t * t * t); }
double ans1(double t) { return (t * t * t * t / 4); }
double func2(double y, double t) { return(-y); }
double ans2(double t) { return (std::sin(t)); }
Eigen::Vector<double, 4> func_orb(const Eigen::Vector<double, 4> state) {
    const double mu = 0.012277471;
    const double eta = 1 - mu;
    double A = std::pow((state[0] + mu) * (state[0] + mu) + state[2] * state[2], 3.0 / 2.0);
    double B = std::pow((state[0] - eta) * (state[0] - eta) + state[2] * state[2], 3.0 / 2.0);
    Eigen::Vector<double, 4> res = { state[1],
                                     state[0] + 2 * state[3] - eta * (state[0] + mu) / A - mu * (state[0] - eta) / B,
                                     state[3],
                                     state[2] - 2 * state[1] - eta * state[2] / A - mu * state[2] / B };
    return (res);
}

СauchyProblem<decltype(func1), 1> linear;
СauchyProblem<decltype(func2), 2> oscillator;
СauchyProblem<decltype(func_orb), 4> ArenstorfOrbit;

int main() {
    /*const double start = 0;
    const double end = 5;

    std::ofstream data1("data_I.txt");
    data1.precision(16);
    std::ofstream data2("data_II.txt");
    data2.precision(16);

    for (double h = -7; h < 0; h += 0.1) {
        
        double step = std::pow(10, h);
        
        std::vector<std::array<double, 2>> res1
            = integrate<decltype(func1), RK4Table, СauchyProblem<decltype(func1), 1>>(func1, Eigen::Vector<double, 1>(0.0), start, end, step, linear);
        double err1 = std::abs(res1[0][0] - ans1(res1[0][1]));

        double delta;
        for (unsigned int i = 1; i < res1.size(); i++) {

            delta = std::abs(res1[i][0] - ans1(res1[i][1]));
            err1 = delta > err1 ? delta : err1;
        }

        std::vector<std::array<double, 2>> res2
            = integrate<decltype(func2), RK4Table, СauchyProblem<decltype(func2), 2>>(func2, Eigen::Vector<double, 2>(0.0, 1.0), start, end, step, oscillator);
        double err2 = std::abs(res2[0][0] - ans2(res2[0][1]));

        for (unsigned int i = 1; i < res2.size(); i++){
            delta = std::abs(res2[i][0] - ans2(res2[i][1]));
            err2 = delta > err2 ? delta : err2;
        }
        
        data1 << h << "     " << err1 << std::endl;
        data2 << h << "     " << err2 << std::endl;
    }

    data1.close();
    data2.close();*/

    std::ofstream data_orb("data_orb.txt");
    data_orb.precision(16);

    const double T = 17.0652165601579625588917206249;

    auto res = integrateDP<decltype(func_orb), DP45, СauchyProblem<decltype(func_orb), 4>>(
        func_orb, { {0.994, 0.0, 0.0, -2.00158510637908252240537862224}, 0 }, 10.0 * T, { 1e-16, 1e-5, 1e-15, 6.3*1e-4 }, ArenstorfOrbit);

    for (int i = 0; i < res.size(); i++) {
        data_orb << res[i].state[0] << "    " << res[i].state[2] << std::endl;
    }

    data_orb.close();
    return 0;
};
