#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <initializer_list>

#include "MatrixLib.h"

using namespace std;
using namespace dml;

template<int cols,int rows>
mat<1, rows, float> GaussianElimination(mat<cols, rows, float> A, mat<1, rows, float> b) {
	upperTriangular(A, b);
	return backSolve(A, b);
}

int main() { 
		typedef mat<4, 4, float> mat4;
		typedef mat<1, 1, float> scalor;
		typedef mat<1, 4, float> vec4;

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

	getchar();
}

