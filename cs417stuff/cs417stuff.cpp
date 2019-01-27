#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <initializer_list>

#include "MatrixLib.h"

using namespace std;
using namespace dml;

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

		upperTriangular(A, b);

		cout << A << endl;
		cout << b << endl;

		vec4 x(0);
		x.data[0][x.numrows() - 1] = b[0][b.numrows() - 1];
		for (int r = 0; r <= A.numrows() - 2; r++) {
			float sum = 0;
			for (int c = r + 1; c < A.numcols() - 1; c++) {
				sum += A[c][r] * b[0][c];
			}
			x.data[0][r] = b[0][r] - sum;
		}
		cout << x << endl;

	getchar();
}

