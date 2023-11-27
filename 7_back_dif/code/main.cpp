#include "FDN.h"
#include<fstream>

Eigen::Vector<double, 1> lin(const Eigen::Vector<double, 1> state, const double t) { 
    Eigen::Vector<double, 1> res;
    res = Eigen::Vector < double, 1>::Zero();
    res[0] = t * t * t;
    return(res);
}
double ans1(double t) { return (t * t * t * t / 4); }
Eigen::Vector<double, 2> osc(const Eigen::Vector<double, 2> state, const double t) {
    Eigen::Vector<double, 2> res = { state[1], -state[0] };
    return(res);
}
double ans2(double t) { return (std::sin(t)); }

Eigen::Vector<double, 6> orbit(const Eigen::Vector<double, 6> state, const double t) {
    Eigen::Vector<double, 6> res;
    const double mu = 398600.441589;
    const double e = 3.0 * mu * 6371.0 * 6371.0 * 1.08263 / 2.0;
    const double r = std::sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
    const double A = -mu / std::pow(r, 3);
    const double B = 5.0 * e / std::pow(r, 7);
    const double C = e / std::pow(r, 5);

    res = { state[3], state[4], state[5],
            A * state[0] + B * state[0] * state[2] * state[2] + C * state[0],
            A * state[1] + B * state[1] * state[2] * state[2] + C * state[1],
            A * state[2] + B * state[2] * state[2] * state[2] - C * state[2] };

    return res;
}

ÑauchyProblem<decltype(lin), 1> linear;
ÑauchyProblem<decltype(osc), 2> oscillator;
ÑauchyProblem<decltype(orbit), 6> Earth;

int main() {
    
    /*const double start = 0;
    const double end = 5;

    std::ofstream data1("data_I.txt");
    data1.precision(16);
    std::ofstream data2("data_II.txt");
    data2.precision(16);

    for (double h = -6; h < 0; h += 0.5) {

        double step = std::pow(10, h);

        auto res1
            = integrate<decltype(lin), BDF4, ÑauchyProblem<decltype(lin), 1>, RK4Table>(lin, { Eigen::Vector<double, 1>(0.0), start }, end, {step, 1e-6, 10}, linear);
        double err1 = (res1[0].state[0] - ans1(start));

        double delta;
        for (unsigned int i = 1; i < res1.size(); i++) {
            
            delta = std::abs(res1[i].state[0] - ans1(res1[i].arg));
            err1 = delta > err1 ? delta : err1;
        }

        auto res2
            = integrate<decltype(osc), BDF4, ÑauchyProblem<decltype(osc), 2>, RK4Table>(osc, { Eigen::Vector<double, 2>(0.0, 1.0), start }, end, { step, 1e-6, 10 }, oscillator);
        double err2 = (res2[0].state[0] - ans2(start));

        for (unsigned int i = 1; i < res2.size(); i++){
            delta = std::abs(res2[i].state[0] - ans2(res2[i].arg));
            err2 = delta > err2 ? delta : err2;
        }

        data1 << h << "     " << err1 << std::endl;
        data2 << h << "     " << err2 << std::endl;
    }

    data1.close();
    data2.close();*/

    std::ofstream data3("data_III.txt");
    data3.precision(16);

    auto res
        = integrate<decltype(orbit), BDF4, ÑauchyProblem<decltype(orbit), 6>, RK4Table>(orbit, { Eigen::Vector<double, 6>(35786, 0.0, 0.0, 0.0, 2.90, 1.0), 0 }, 100000, { 1, 1e-4, 20 }, Earth);

    for (unsigned int i = 0; i < res.size(); i++) {
        data3 << res[i].state[0] << "   " << res[i].state[1] << "   " << res[i].state[2] << std::endl;
    }

    data3.close();

    return 0;
};
