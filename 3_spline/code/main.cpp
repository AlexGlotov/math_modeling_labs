#include"funcs.h"
#include<fstream>

int main() {
    
    const double start = 0;
    const double end = 10;
    
    std::ofstream data1("data1.txt");
    data1.precision(10);

    for (unsigned int N = 5; N < 1000; N++) {
        std::vector<double> points = linear_grid<double>(start, end, N);
        std::vector<double> values;

        for (unsigned int i = 0; i < N; i++) {
            values.push_back(std::exp(points[i]));
        }

        CubicSpline<double, double> spline(points, values);

        double err = 0;
        for (double x = start; x < end; x += end / 10000) {

            double delta = abs(spline.interpolate(x) - exp(x));
            err = delta > err ? delta : err;
        }

        data1 << N << "      " << err << std::endl;
    }
    
    data1.close();
    return 0;
}