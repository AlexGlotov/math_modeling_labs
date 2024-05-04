#include"funcs.h"
#include<fstream>
#include<math.h>

const double pi = acos(-1.0);
double f(double x, double t)
{
    return t * (x + 1);
}

double start_cond(double x)
{
    return 0;
}

double F_left(double t)
{
    return t * t;
}

double F_right(double t)
{
    return t * t;
}



int main() {

    const double time_start = 0;
    const double time_end = 10;
    const double tau = 0.01;

    const double x_start = 0;
    const double x_end = 1.0;

    BoundCond consts{};
    consts.a_left = 0;
    consts.a_right = 1;
    consts.b_left = 1;
    consts.b_right = 0;

    std::ofstream data("data7.txt");
    data.precision(10);

    const unsigned int N = 6999;
    double step = (x_end - x_start) / N;
    const double CFL = 1.0 * tau / (2 * step * step);
    data << tau << "    " << step << "      " << N << std::endl;

    std::vector<double> res = Create_Mesh(x_start, x_end, step, N, start_cond);

    for (double t = time_start; t < time_end; t += tau) {

        data << t << "		";
        for (int i = 0; i < res.size(); i++)
        {
            data << res[i] << "		";
        }
        data << std::endl;

        std::vector<double> column = create_column(res, CFL, step, f, F_left, F_right, x_start,
               tau, t, consts.b_left, consts.b_right);

        ThreeDiagonalMatrix<double> matrix(step, consts, CFL, N);

        res = solve(matrix, column);
        
    }

    return 0;
}