#include "Header.h"
#include<fstream>

cell u_0(double x)
{
	if (x <= 0)
	{
		cell c;
		c.x = x;
		c.rho = 13.0;
		c.u = 0.0;
		c.e = 1500000.0 / 13.0;
		c.p = 1000000.0;
		c.w = Eigen::Vector<double, 3>{ 13.0, 0.0 , 1500000.0 };
		return (c);
	}
	else
	{
		cell c;
		c.x = x;
		c.rho = 1.3;
		c.u = 0.0;
		c.e = 150000.0 / 1.3;
		c.p = 100000.0;
		c.w = Eigen::Vector<double, 3>{ 1.3, 0.0 , 150000.0 };
		return (c);
	};
}

int main()
{
	const double T = 0.02;
	const double h = 0.1;
	const double L = 10.0;
	const double tau = 1e-5;

	std::ofstream data_rho("data_rho.txt");
	std::ofstream data_u("data_u.txt");
	std::ofstream data_e("data_e.txt");
	std::ofstream data_p("data_p.txt");
	data_rho.precision(16);
	data_u.precision(16);
	data_e.precision(16);
	data_p.precision(16);

	mesh1D res = Create_Mesh(L, h, u_0);

	double t = 0;

	data_rho << "rho-" << "	";
	data_u << "u-" << "	";
	data_e << "e-" << "	";
	data_p << "p-" << "	";
	for (int i = 0; i < res.layer.size(); i++) {
		data_rho << res.layer[i].x << "		";
		data_u << res.layer[i].x << "		";
		data_e << res.layer[i].x << "		";
		data_p << res.layer[i].x << "		";
	}
	data_rho << std::endl;
	data_u << std::endl;
	data_e << std::endl;
	data_p << std::endl;

	while (t < T){

		t = res.t;
		data_rho << res.t << "	";
		data_u << res.t << "	";
		data_e << res.t << "	";
		data_p << res.t << "	";
		for (int i = 0; i < res.layer.size(); i++)
		{
			data_rho << res.layer[i].rho << "		";
			data_u << res.layer[i].u << "		";
			data_e << res.layer[i].e << "		";
			data_p << res.layer[i].p << "		";
		}
		data_rho << std::endl;
		data_u << std::endl;
		data_e << std::endl;
		data_p << std::endl;

		res = KIR_scheme(res, tau, h, 5.0 / 3.0);
	}

	return 0;
}