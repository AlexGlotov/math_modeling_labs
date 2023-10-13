#include"funcs.h"
#include<fstream>

int main() {
    
    const double start = 0;
    const double end = 10;
    
    std::ofstream data1("data2.txt");
    data1.precision(10);

    for (unsigned int N = 5; N < 1000; N++) {
        std::vector<double> points = linear_grid<double>(start, end, N);
        std::vector<double> values;

        for (unsigned int i = 0; i < N; i++) {
            values.push_back(std::exp(points[i]));
        }

        double left_der = std::exp(start), right_der = std::exp(end); 
     
     
        CubicSpline<double, double> spline(points, values, left_der, right_der);

        double err = 0;
        for (double x = start; x < 0.999 * end; x += end / 10000) {

            double delta = abs(spline.interpolate(x) - exp(x));
            err = delta > err ? delta : err;
        }

        data1 << N << "      " << err << std::endl;
    }

    data1.close();

    return 0;
}