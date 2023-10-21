#include "funcs.h"
#include<fstream>
#include<string>

template <typename Type> double sgn(Type x) {
	return (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0);
}

double func_x(double x) { return (x * x + std::tan(x) * std::tan(x) - 1); }
double func_y(double x) { return (x - sgn(x) * std::tan(std::sqrt(1 - x * x))); }

int main() {
	/*const unsigned int N = 4;

	const double M = std::acos(-1) / 4;
	std::array<double, N>e{ 0.1, 0.2, 0.5, 0.8 };
	std::array<double, N>ans{ 0.861264884868, 0.9478282237996, 1.261703055253, 1.585313861354 };

	for (int i = 0; i < N; i++) {
		std::string name = "data" + std::to_string(i) + ".txt";
		std::ofstream data(name);
		data.precision(16);

		for (int j = 1; j < 10; j++) {
			data << j << "		" << std::log10(std::abs(ans[i] -
				keplerSolver(e[i], M, j, 0))) << std::endl;
		}

		data.close();
	}*/
	std::cout.precision(7);
	std::cout << "(" << solve(func_x, -0.1861978, 0.5, 15) <<
		"	" << solve(func_y, -0.339661, 0.5, 10) << ")" << std::endl;
	std::cout << "(" << solve(func_x, 0.1861978, -0.5, 15) <<
		"	" << solve(func_y, -0.339661, -0.5, 10) << ")" << std::endl;

	return 0;
}