#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <initializer_list>

#include "MatrixLib.h"

using namespace std;
using namespace dml;

template<int cols,int rows,class TYPE>
mat<1, rows, TYPE> GaussianElimination(mat<cols, rows, TYPE> A, mat<1, rows, TYPE> b) {
	upperTriangular(A, b);
	return backSolve(A, b);
}

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
			{3},
			{2},
			{9},
			{4}});

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



	getchar();
}

