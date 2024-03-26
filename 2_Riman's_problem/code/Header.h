#pragma once

#include<vector>
#include<iostream>
#include<math.h>
#include"Eigen/Dense"

struct cell {
    double x;
    double rho; //density
    double u; //velocity
    double e; //energy
    double p; //pressure
    Eigen::Vector<double, 3> w = { rho, rho * u, rho * e };
};

struct mesh1D {
    std::vector<cell> layer;
    double t;
};

template<typename Callable>
mesh1D Create_Mesh(const double L, const double h, const Callable& u_0) {
    std::vector<cell> data;
    const unsigned int N = L / h + 1;
    const double step = L / N;

    cell c;
    c = u_0(-L);
    data.push_back(c);
    for (unsigned int i = 1; i < 2 * N + 1; i++) {
        c = u_0(data[i - 1].x + step);
        data.push_back(c);
    }

    return { data, 0.0 };
}

double max_abs_lambda(const double u, const double c) {
    double M = (abs(c - u) > c) ? abs(c - u) : c;
    M = (M > abs(c + u)) ? M : abs(c + u);
    return M;
}

/*the transfer equation, the Kurrent-Izakson-Rees scheme*/
mesh1D KIR_scheme(mesh1D input_data, double tau, const double h, const double gamma)
{
    std::vector<cell> data = input_data.layer;
    double t = input_data.t;
    unsigned int N = (input_data.layer).size();
    
    double max_lambda = 0;
    for (int i = 0; i < N; i++) {
        double c = sqrt(gamma * (gamma - 1) * data[i].e);
        max_lambda = (max_lambda > max_abs_lambda(data[i].u, c)) 
            ? max_lambda : max_abs_lambda(data[i].u, c);
    }

    while ((tau * max_lambda / h) > 0.01) {
        tau = tau / 2.0;
    }
    t += tau;

    for (int i = 2; i < data.size() - 2; i++) {
        double u = data[i].u;
        double c = sqrt(gamma * (gamma - 1) * data[i].e);

        Eigen::Matrix < double, 3, 3 > Omega{
                                       {-u * c, c, gamma - 1},
                                       {-c * c, 0, gamma - 1},
                                       {u * c, -c, gamma - 1} };

        Eigen::Matrix < double, 3, 3 > Omega_rev{
                                           {1 / (2 * c * c), -1 / (c * c), 1 / (2 * c * c)},
                                           {(u + c) / (2 * c * c), -u / (c * c),(u - c) / (2 * c * c)},
                                           {1 / (2 * gamma - 2), 0, 1 / (2 * gamma - 2)} };

        Eigen::Matrix<double, 3, 3> Lambda{
                                          {u + c, 0 , 0},
                                          {0, u, 0},
                                          {0, 0, u - c} };

        Eigen::Matrix<double, 3, 3> abs_Lambda{
                                              {abs(u + c), 0, 0},
                                              {0, abs(u), 0},
                                              {0, 0, abs(u - c)} };

        data[i].w = data[i].w - tau * Omega_rev * Lambda * Omega / (2 * h) * (data[i + 1].w - data[i - 1].w) +
            tau * Omega_rev * abs_Lambda * Omega / (2 * h) * (data[i + 1].w - 2 * data[i].w + data[i - 1].w);

        data[i].rho = data[i].w[0];
        data[i].u = data[i].w[1] / data[i].rho;
        data[i].e = data[i].w[2] / data[i].rho;
        data[i].p = (gamma - 1) * data[i].rho * data[i].e;
    }

    

    return { data, t };
}
