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
		typedef mat<4, 4, double> mat4;
		typedef mat<1, 4, double> vec4;

		mat4 A=singularMatrix<4,4,double>(1,10,1,true);
		vec4 b=randomVector<4,double>(1,10,true);

		auto x = GaussianElimination(A, b);
		auto z = A * x;
		cout << "b:\n";
		cout << b << endl;
		cout << "z:\n";
		cout << z << endl;

		mat4 L, U;
		decomposeMatrix(A, L, U);

		cout << "L:\n";
		cout << L << endl;
		cout << "U:\n";
		cout << U << endl;

		auto out = U*L;
		cout << "A:\n";
		cout << A << endl;
		cout << "LU:\n";
		cout << out << endl;
		cout << endl<<endl;

	getchar();
}

