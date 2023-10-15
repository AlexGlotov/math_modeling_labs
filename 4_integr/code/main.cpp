#include "funcs.h"
#include <fstream>

double func(double x) { return std::sin(x); }

int main() 
{
	const double start = 0;
	const double end = 10;

	std::ofstream data1("data1.txt");
	data1.precision(16);
	std::ofstream data2("data2.txt");
	data2.precision(16);
	std::ofstream data3("data3.txt");
	data2.precision(16);
	std::ofstream data4("data4.txt");
	data2.precision(16);
	std::ofstream data5("data5.txt");
	data2.precision(16);

	const double I = 1 - std::cos(10);
	for (Dif<typename ArgumentGetter<double(double)>::Argument> h = -7; h < 1; h += 0.001) {
		double step = std::exp(h);
		data1 << h << "		" << std::abs(I - integrate<decltype(func), double, 1>
			(func, start, end, nodes<double, 1>::p, nodes<double, 1>::w, step)) << std::endl;

		data2 << h << "		" << std::abs(I - integrate<decltype(func), double, 2>
			(func, start, end, nodes<double, 2>::p, nodes<double, 2>::w, step)) << std::endl;

		data3 << h << "		" << std::abs(I - integrate<decltype(func), double, 3>
			(func, start, end, nodes<double, 3>::p, nodes<double, 3>::w, step)) << std::endl;

		data4 << h << "		" << std::abs(I - integrate<decltype(func), double, 4>
			(func, start, end, nodes<double, 4>::p, nodes<double, 4>::w, step)) << std::endl;

		data5 << h << "		" << std::abs(I - integrate<decltype(func), double, 5>
			(func, start, end, nodes<double, 5>::p, nodes<double, 5>::w, step)) << std::endl;

	}

	data1.close();
	data2.close();
	data3.close();
	data4.close();
	data5.close();

	return 0;
}