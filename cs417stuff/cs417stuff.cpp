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

int menu() {
	int option;
	cout << "  Menu\n";

	cin >> option;
	return option;
}
void interactiveMenu() {

	int option = -1;
	while (option != 0) {
		option = menu();
		switch (option){
		case 1: {}break;

		}
	}
}

int main() {
	const size_t N = 4;
	double min = 0, max = 100;
	default_random_engine gen(time(NULL));

	mat A = nonsingularMatrix(gen,min, max, N, true);
	mat b = randomVector(gen, min, max,N, true);

	auto x = GaussianElimination(A, b);
	auto z = A * x;
	cout << "b:\n";
	cout << b << endl;
	cout << "Ax:\n";
	cout << z << endl;

	mat L, U;
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
		mat At = nonsingularMatrix(gen, min, max, N, true);
		mat bt = randomVector(gen, min, max,N, true);
		vector<double> error;
		auto xnew = GaussSeidel(At, bt, randomVector(gen, min, max,N, true),error,1, 1000);
		auto twoNorm = At * xnew - bt;
		if(norm(twoNorm)==0)
			cout << error.size() << endl<<endl;
	}

	getchar();
}

