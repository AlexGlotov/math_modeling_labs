#include "Header.h"
#include<fstream>

double u_0(double x)
{
	if ((x >= 0.0) and (x <= 5.0)) {
		return sin(x * 3.1415 / 5.0);
	}
	else return 0.0;
}

int main()
{
	const double T = 18;
	const double h = 0.5;
	const double CFL = 1.01;
	const double L = 20.0;
	const double tau = CFL * h;
	const double a = 1.0;

	std::ofstream data("data-LW-1.01.txt");
	data.precision(16);

	data << T << "	" << h << std::endl;

	mesh1D res = Create_Mesh(L, h, u_0);

	for (double t = 0; t < T; t += tau) {
		
		data << res.t << "		";
		for (int i = 0; i < res.layer.size(); i++)
		{
			data << res.layer[i].u << "		";
		}
		data << std::endl;

		//res = KIR_scheme(res, a, tau, h);
		res = LW_scheme(res, a, tau, h);
	}

	return 0;
}