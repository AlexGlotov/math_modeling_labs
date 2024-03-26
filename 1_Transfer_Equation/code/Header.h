#pragma once

#include<vector>
#include<iostream>
#include<math.h>


struct cell {
    double x;
    double u;
};

struct mesh1D {
    std::vector<cell> layer;
    double t;
};

template<typename Callable>
mesh1D Create_Mesh(const double L, const double h, const Callable& u_0) {
    std::vector<cell> data;
    const unsigned int N = 20;

    cell c;
    c.x = 0;
    c.u = u_0(0);
    data.push_back(c);
    for (unsigned int i = 1; i < N; i++) {
        cell c;
        c.x = data[i - 1].x + h;
        c.u = u_0(c.x);
        data.push_back(c);
    }

    return { data, 0.0 };
}


/*the transfer equation, the Kurrent-Izakson-Rees scheme*/
mesh1D KIR_scheme(mesh1D input_data, const double a, const double tau, const double h)
{
    mesh1D data = input_data;
    unsigned int N = (input_data.layer).size();
    data.t += tau;
    for (unsigned int i = 1; i < (data.layer).size() - 1; i++) {
        data.layer[i].u += (tau / (2 * h)) * (input_data.layer[i - 1].u * (a + abs(a)) -
            2 * abs(a) * input_data.layer[i].u + (abs(a) - a) * input_data.layer[i + 1].u);
    }

    data.layer[0].u += (tau / (2 * h)) * ((input_data).layer[N - 1].u * (a + abs(a)) -
        2 * abs(a) * (input_data).layer[0].u + (abs(a) - a) * (input_data).layer[1].u);

    data.layer[N - 1].u += (tau / (2 * h)) * ((input_data).layer[N - 2].u * (a + abs(a)) -
        2 * abs(a) * (input_data).layer[N - 1].u + (abs(a) - a) * (input_data).layer[0].u);

    return data;
}

/*the transfer equation, the Lax-Wendroff scheme*/
mesh1D LW_scheme(mesh1D input_data, const double a, const double tau, const double h)
{
    mesh1D data = input_data;
    unsigned int N = (input_data.layer).size();
    data.t += tau;
    const double Co = a * tau / h;
    for (unsigned int i = 1; i < (data.layer).size() - 1; i++) {
        data.layer[i].u = (Co * Co + Co) * input_data.layer[i - 1].u / 2.0 +
            (1 - Co * Co) * input_data.layer[i].u + (Co * Co - Co) * input_data.layer[i + 1].u / 2.0;
    }

    data.layer[0].u = (Co * Co + Co) * input_data.layer[N - 1].u / 2.0 +
        (1 - Co * Co) * input_data.layer[0].u + (Co * Co - Co) * input_data.layer[1].u / 2.0;

    data.layer[N - 1].u = (Co * Co + Co) * input_data.layer[N - 2].u / 2.0 +
        (1 - Co * Co) * input_data.layer[N - 1].u + (Co * Co - Co) * input_data.layer[0].u / 2.0;

    return data;
}