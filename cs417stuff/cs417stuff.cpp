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

		mat4 A({ 
			{1,1,3,7}, 
			{4,2,9,4}, 
			{3,4,3,4}, 
			{6,2,5,4}, 
		});
		vec4 b({
			{1},
			{3},
			{1},
			{9}});

		auto x = GaussianElimination(A, b);
		auto z = A * x;
		cout << z << endl;
		cout << b << endl;

		mat4 L, U;
		decomposeMatrix(A, L, U);

		cout << L << endl;
		cout << U << endl;

		auto out = U*L;
		cout << out << endl;
		cout << A << endl;
		cout << endl<<endl;
		mat4 m = singularMatrix<4,4,double>(1);
		cout << m << endl;
		vec4 bb({
			{ 6 },
			{ 8 },
			{ 2 },
			{ 7 },
			});
		auto xx = GaussianElimination(m, bb);
		auto zz = m * xx;
		cout << bb << endl;
		cout << zz << endl;

	getchar();
}

