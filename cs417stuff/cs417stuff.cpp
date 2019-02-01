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
	const size_t N = 4;
	double min = 0, max = 100;
	default_random_engine gen(time(NULL));

	typedef mat<N, N, double> mat4;
	typedef mat<1, N, double> vec4;

	mat4 A = nonsingularMatrix<N, N, double>(gen,min, max, 1, true);
	vec4 b = randomVector<N, double>(gen, min, max, true);

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

	for (int i = 0; i < 100; i++) {
		mat4 At = nonsingularMatrix<N, N, double>(gen, min, max, 1, true);
		vec4 bt = randomVector<N, double>(gen, min, max, true);
		vector<double> error;
		auto xnew = JacobiIterative(At, bt, randomVector<N, double>(gen, min, max, true),error,1, 1000);
		auto twoNorm = At * xnew - bt;
		if(norm(twoNorm)==0)
			cout << error.size() << endl<<endl;
	}

	getchar();
}

