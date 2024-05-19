#include<fstream>
#include<iostream>
#include<cmath>
#include"Eigen/Dense"
#include<vector>

const double PI = acos(-1);
const double lambda = 0.0001;

double u0t(double x, double y, double t)
{
	return (cos(PI * x) * sin(5 * PI * y));
}

double u0y(double x, double y, double t)
{
	return 0;
}

double u1y(double x, double y, double t)
{
	return 0;
}

double u0x(double x, double y, double t)
{
	return (sin(5 * PI * y) * exp(-50 * PI * PI * lambda * t));
}

double u1x(double x, double y, double t)
{
	return (-sin(5 * PI * y) * exp(-50 * PI * PI * lambda * t));
}

double ans(double x, double y, double t)
{
	return (cos(PI * x) * sin(5 * PI * y) * exp(-50 * PI * PI * lambda * t));
}

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

class ThreeDiagonalMatrix {
private:
	std::vector<std::array<double, 3>> diag;
public:
	ThreeDiagonalMatrix(const double CFL, const int N) {

		diag.push_back({ 0.0, 1.0 , 0.0 });
		for (unsigned int i = 1; i < N; i++)
		{
			diag.push_back({ -CFL, 1.0 + 2.0 * CFL, -CFL });
		}
		diag.push_back({ 0.0, 1.0, 0.0 });

	}

	std::array<double, 3> operator() (const unsigned int i) const {
		return diag[i];
	};

};

std::vector<double> create_column(std::vector<double> data, const double CFL, double left_value, double right_value)
{
	std::vector<double> col;
	unsigned int S = data.size();

	col.push_back(0);
	for (int i = 1; i < S - 1; i++)
	{
		col.push_back(data[i]);

	}

	col.push_back(right_value);

	col[0] = left_value;

	return col;
}

template<unsigned int N>
std::vector<double> solve(const ThreeDiagonalMatrix& matrix,
	const std::vector<double>& column)
{
	unsigned int N = column.size();
	int n = N - 1;
	std::vector<double> p_vector{ 0 }, q_vector{ 0 };
	p_vector.reserve(N);
	q_vector.reserve(N);
	for (int i = 0; i < N; ++i) {
		p_vector.push_back(-matrix(i)[2] /
			(matrix(i)[0] * p_vector[i] + matrix(i)[1]));
		q_vector.push_back((column[i] - matrix(i)[0] * q_vector[i]) /
			(matrix(i)[0] * p_vector[i] + matrix(i)[1]));
	}
	std::vector<double> sol(N);
	sol[n] = (column[n] - matrix(n)[0] * q_vector[n]) /
		(matrix(n)[0] * p_vector[n] + matrix(n)[1]);
	for (int i = n - 1; i >= 0; --i) {
		sol[i] = p_vector[i + 1] * sol[i + 1] + q_vector[i + 1];
	}
	return sol;
};


template<typename Callable, const unsigned int Nx, const unsigned int Ny>
double calc(const unsigned int N, Callable u0t, Callable u0x, Callable u1x, Callable u0y, Callable u1y, Callable ans,
	const double a_x, const double a_y)
{
	std::vector<double> tt;
	const double tau = 0.5 / N;
	tt.push_back(0.0);
	std::vector<double> xx;
	const double hx = 1.0 / Nx;
	xx.push_back(0.0);
	std::vector<double> yy;
	const double hy = 1.0 / Ny;
	yy.push_back(0.0);
	for (unsigned int i = 1; i <= 2 * N; i++) tt.push_back(tt[i - 1] + tau);
	for (unsigned int i = 1; i <= Nx; i++) xx.push_back(xx[i - 1] + hx);
	for (unsigned int i = 1; i <= Ny; i++) yy.push_back(yy[i - 1] + hy);

	Eigen::MatrixXd u(Nx + 1, Ny + 1);
	Eigen::MatrixXd u_1_2(Nx + 1, Ny + 1);

	for (unsigned int i = 0; i <= Nx; i++) {
		for (unsigned int j = 0; j <= Ny; j++) {
			u(j, i) = u0t(xx[i], yy[j], 0);
		}
	}

	double CoY = a_y * 2 * tau * lambda / (hy * hy);
	double CoX = a_x * 2 * tau * lambda / (hx * hx);

	ThreeDiagonalMatrix Ax(CoX, Nx);
	ThreeDiagonalMatrix Ay(CoY, Ny);

	for (unsigned int n = 0; n < N; n++) {

		for (unsigned int i = 0; i <= Ny; i++) {
			std::vector<double> data;
			for (unsigned int j = 0; j <= Nx; j++) data.push_back(u(i, j));
			std::vector<double> d = create_column(data, CoX, u0x(0, yy[i], tt[2 * n + 1]), u1x(1, yy[i], tt[2 * n + 1]));
			std::vector<double> sol = solve<Nx + 1>(Ax, d);
			for (unsigned int j = 0; j <= Nx; j++) u_1_2(i, j) = sol[j];
		}
		for (unsigned int i = 0; i <= Nx; i++) {
			std::vector<double> data;
			for (unsigned int j = 0; j <= Ny; j++) data.push_back(u_1_2(j, i));
			std::vector<double> d = create_column(data, CoY, u0y(xx[i], 0, tt[2 * (n + 1)]), u1y(xx[i], 1, tt[2 * (n + 1)]));
			std::vector<double> sol = solve<Ny + 1>(Ay, d);
			for (unsigned int j = 0; j <= Ny; j++) u(j, i) = sol[j];
		}
	}

	Eigen::MatrixXd answ(Nx + 1, Ny + 1);
	Eigen::MatrixXd delta(Nx + 1, Ny + 1);
	for (unsigned int i = 0; i <= Nx; i++) {
		for (unsigned int j = 0; j <= Ny; j++) {
			answ(j, i) = ans(xx[i], yy[j], 1.0);
		}
	}
	delta = u - answ;
	double err = abs(delta.maxCoeff()) > abs(delta.minCoeff()) ? abs(delta.maxCoeff()) : abs(delta.minCoeff());
	return err;
}

int main()
{
	const unsigned int Nx = 50;
	const unsigned int Ny = 50;
	unsigned int N = 1000;

	std::cout << "Nt = " << N << " , Nx = " << Nx << " , Ny = " << Ny << " , err = " << calc<decltype(u0t), Nx, Ny>(N, u0t, u0x, u1x, u0y, u1y, ans, 25.0, 1.0);


	return 0;
}