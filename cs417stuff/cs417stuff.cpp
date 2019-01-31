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
	const size_t N = 8;
	double min = 0, max = 100;

	typedef mat<N, N, double> mat4;
	typedef mat<1, N, double> vec4;

	mat4 A = singularMatrix<N, N, double>(min, max, 1, true);
	vec4 b = randomVector<N, double>(min, max, true);

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

	vector<double> error;
	auto xnew = JacobiIterative(A, b, randomVector<N, double>(min, max, true),error,1, 100);
	auto twoNorm = A * xnew - b;
	cout << x << endl << xnew << endl;
	cout << norm(twoNorm) << endl<<endl;

	for (int i = 0; i < error.size(); i++) {
		cout << error[i] << endl;
	}

	getchar();
}

