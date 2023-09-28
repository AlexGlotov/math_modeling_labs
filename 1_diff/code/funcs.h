#pragma once

#include<iostream>
#include<array>
#include<eigen-3.4.0\Eigen\Dense>

using namespace std;

int factorial(int i)
{
	if (i == 0) return 1;
	else return i * factorial(i - 1);
};

template<typename RealType, unsigned int N>
struct DerivativeCoef {
	RealType centralCoef;
	array<RealType, N> otherCoeffs;
};

template<typename RealType, unsigned int N, unsigned int L>
DerivativeCoef<RealType, N>
calcDerivativeCoef(const array<RealType, N>& points) noexcept {
	Eigen::Matrix<RealType, N + 1, N + 1> A = Eigen::Matrix<RealType, N + 1, N + 1>::Ones();
	for (unsigned int i = 1; i <= N; i++) {
		A(i, 0) = 0;
		for (unsigned int j = 1; j <= N; j++) {
			A(i, j) = A(i - 1, j) * points[j - 1];
		}
	}

	Eigen::Matrix<RealType, N + 1, 1> b = Eigen::Matrix<RealType, N + 1, 1>::Zero();
	b(L) = factorial(L);

	Eigen::Vector<RealType, N + 1> c = A.colPivHouseholderQr().solve(b);

	RealType centralCoef = c(0);

	array<RealType, N> coefs;
	for (unsigned int i = 0; i < N; i++) {
		coefs[i] = c(i + 1);
	}

	return DerivativeCoef<RealType, N>{ centralCoef, coefs };
}

template<typename RealType, unsigned int N, unsigned int L>
RealType d_e(const RealType x0, const RealType h, const array<RealType, N>& points) {
	DerivativeCoef<double, N> coefs = calcDerivativeCoef<RealType, N, L>(points);
	double order_const = pow(h, L);

	RealType de_x0 = coefs.centralCoef * exp(x0) / order_const;
	for (unsigned int i = 0; i < N; i++) {
		de_x0 += coefs.otherCoeffs[i] * exp(x0 + points[i] * h) / order_const;
	}

	return de_x0;
}
