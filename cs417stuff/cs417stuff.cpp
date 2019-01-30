#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <initializer_list>
#include <random>
#include <ctime>

#include "MatrixLib.h"

using namespace std;
using namespace dml;



int main() {
	const int N = 4;
	double min = 0, max = 100;

	typedef mat<N, N, double> mat4;
	typedef mat<1, N, double> vec4;

	mat4 A = singularMatrix<N, N, double>(min, max, 1,true);
	vec4 b = randomVector<N, double>(min, max,true);
	{
		mat4 L, D, U;
		split(A, L, D, U);
		auto out = L + D + U;
		cout << A << endl << out << endl;
		cout << endl << endl;

		mat4 Dinv = 1. / D;//make its column vec
		mat4 LU = L + U;
		vec4 xold= randomVector<N, double>(min, max, true), xnew;
		int iters = 1;
		for (int i = 0; i < iters; i++) {
			xnew = Dinv * (b - LU * xold);
			xold = xnew;
		}
		cout << b << endl << xnew << endl;
	}
	/*{
		auto x = GaussianElimination(A, b);
		auto z = A * x;
		cout << "b:\n";
		cout << b << endl;
		cout << "Ax:\n";
		cout << z << endl;

		mat4 L, U;
		decomposeMatrix(A, L, U);

		cout << "L:\n";
		cout << L << endl;
		cout << "U:\n";
		cout << U << endl;

		auto out = L * U;
		cout << "A:\n";
		cout << A << endl;
		cout << "LU:\n";
		cout << out << endl;

		auto res = U * x;
		auto y = forwordSolve(L, b);
		cout << "y:\n" << y << endl;
		cout << "Ux:\n" << res << endl;

	}*/
	getchar();
}

