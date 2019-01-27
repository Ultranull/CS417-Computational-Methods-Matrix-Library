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

		vec4 x=backSolve(A,b);
		cout << x << endl;

	getchar();
}

