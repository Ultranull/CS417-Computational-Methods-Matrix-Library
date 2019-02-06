#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <chrono>
#include <string>

#include "MatrixLib.h"

using namespace std;
using namespace dml;
using namespace std::chrono;

void load(mat &A, mat &b,string file) {
	ifstream fin;
	fin.open(file.c_str(), ios::in);
	int n;
	fin >> n;
	A = mat(n, n);
	b = mat(1, n);
	double val;
	for(int i=0;i<n;i++)
		for (int j = 0; j < n; j++) {
			fin >> val;
			A[j][i] = val;
	}
	for (int i = 0; i < n; i++) {
		fin >> val;
		b[0][i] = val;
	}
	fin.close();
}

void save(mat &A, mat &b, string file) {
	ofstream fout;
	fout.open(file.c_str(), ios::out);
	int n=A.cols();
	fout << n<<endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << A[j][i] << " ";
		}
		fout << endl;
	}
	for (int i = 0; i < n; i++) {
		fout << b[0][i]<<" ";
	}
	fout.close();
}

int menu() {
	int option;
	cout << "\n    Menu\n";
	cout << "1. generate random matrix A and right side b\n";
	cout << "2. Load matrix A and vector b from a file\n";
	cout << "3. Save current matrix A and vector b to a file.\n";
	cout << "4. Solve current matrix using Gaussian Elimination\n";
	cout << "5. Solve current matrix using LU decomposition\n";
	cout << "6. Solve current matrix with Jacobi's iterative method\n";
	cout << "7. Solve current matrix with Gauss-Seidel method\n";
	//cout << "9. Save current solution vector x to file\n";
	//cout << "10. Graph current solution vector x using GNUPLOT\n";
	cout << "13. display A and b\n";
	cout << "0. exit\n";

	cin >> option;
	return option;
}
void interactiveMenu() {
	default_random_engine gen(time(NULL)); 
	high_resolution_clock::time_point t1,t2;
	mat A, b,L,U,D,x;
	int option = -1;
	while (option != 0) {
		option = menu();
		switch (option){
		case 1: {
			int n;
			cout << "what size?: ";
			cin >> n;
			A = nonsingularMatrix(gen,0,10,n,true);
			b = randomVector(gen,0,10,n,true);
		}break;
		case 2: {
			string n;
			cout << "what file?: ";
			cin >> n;
			load(A, b, n);
		}break;
		case 3: {
			string n;
			cout << "what file?: ";
			cin >> n;
			save(A, b, n);
		}break;
		case 4: {
			t1 = high_resolution_clock::now();
			x = GaussianElimination(A, b);
			cout << "\nx:\n" << x;
			cout << duration<double, std::milli>(high_resolution_clock::now() - t1).count() << " milliseconds \n";
		}break;
		case 5: {
			t1 = high_resolution_clock::now();
			x = LUDecomposition(A, b);
			cout << "\nx:\n" << x;
			cout << duration<double, std::milli>(high_resolution_clock::now() - t1).count() << " milliseconds \n";
		}break;
		case 6: {
			int iters;
			cout << "number of iterations?: ";
			cin >> iters;
			t1 = high_resolution_clock::now();
			x = JacobiIterative(A, b, b, iters);
			cout << "\nx:\n" << x;
			cout << "2-norm: " << norm(A * x - b)<<endl;
			cout << duration<double, std::milli>(high_resolution_clock::now() - t1).count() << " milliseconds \n";
		}break;
		case 7: {
			int iters;
			cout << "number of iterations?: ";
			cin >> iters;
			t1 = high_resolution_clock::now();
			x = GaussSeidel(A, b, b, iters);
			cout << "\nx:\n" << x;
			cout << "2-norm: " << norm(A * x - b) << endl;
			cout << duration<double, std::milli>(high_resolution_clock::now() - t1).count() << " milliseconds \n";
		}break;
		case 13: {
			cout << "\nA: \n" << A << "b:\n" << b << endl;
		}break;
		}


	}
}

int main() {
	{
		//const size_t N = 4;
		//double min = 0, max = 100;
		//default_random_engine gen(time(NULL));
		///*
		//mat A = nonsingularMatrix(gen, min, max, N, true);
		//mat b = randomVector(gen, min, max, N, true);

		//auto x = GaussianElimination(A, b);
		//auto z = A * x;
		//cout << "b:\n";
		//cout << b << endl;
		//cout << "Ax:\n";
		//cout << z << endl;

		//mat L, U;
		//decomposeMatrix(A, L, U);

		//cout << "L:\n";
		//cout << L << endl;
		//cout << "U:\n";
		//cout << U << endl;

		//mat out = L * U;
		//cout << "A:\n";
		//cout << A << endl;
		//cout << "LU:\n";
		//cout << out << endl;

		//mat res = U * x;
		//mat y = forwordSolve(L, b);
		//cout << "y:\n" << y << endl;
		//cout << "Ux:\n" << res << endl;
		//*/
		//for (int i = 0; i < 100; i++) {
		//	mat At = nonsingularMatrix(gen, min, max, N, true);
		//	mat bt = randomVector(gen, min, max, N, true);
		//	vector<double> error;
		//	mat xnew = SOR(At, bt, bt,1.2, error, 1, 100);
		//	mat twoNorm = At * xnew - bt;
		//	cout << "SOR: " <<error.size()<<" 2norm: "<<norm(twoNorm) << endl;
		//	error.clear();
		//	xnew = GaussSeidel(At, bt, bt, error, 1, 100);
		//	twoNorm = At * xnew - bt;
		//	cout <<"guass-siedel: "<< error.size() << " 2norm: " << norm(twoNorm) << endl << endl;
		//}
	}
	interactiveMenu();
	cout << "done!\n";
	getchar();
}

