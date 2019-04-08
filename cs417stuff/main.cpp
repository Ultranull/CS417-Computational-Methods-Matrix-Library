#include <vector>
#include <iostream>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <chrono>
#include <string>
#include <stdlib.h>

#define MAT_IMPL
#include "MatrixLib.h"
#define POLY_IMPL
#include "PolynomLib.h"

using namespace std;
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
	cout << "8. Solve current matrix with SOR method\n";
	cout << "9. Find EiganVector and dominant Eigan Value using power iteration\n";
	cout << "10. Find a root to polynomial using newtons method\n";
	cout << "11. Find a root to polynomial using bisection\n";
	cout << "12. Graph current solution vector x using GNUPLOT\n";
	cout << "13. display A and b\n";
	cout << "14. display stored data information\n";
	cout << "0. exit\n";

	cin >> option;
	return option;
}
void interactiveMenu() {
	default_random_engine gen(time(NULL)); 
	high_resolution_clock::time_point t1;
	mat A(1,2,2), b(1, 1, 2),L,U,D,x,eiganvec;
	double eiganVal=-0xff;
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
			x = mat();
			eiganvec = mat();
			eiganVal = -0xff;
		}break;
		case 2: {
			string n;
			cout << "what file?: ";
			cin >> n;
			load(A, b, n);
			x = mat();
			eiganvec = mat();
			eiganVal = -0xff;
			if (A.rows() >= 300)
				cout << "\n!iterative method is recomended for matrix size!\n(use Gauss-Seidel with -1 (infinite iterations) for best results)\n";
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
			cout << duration<double, std::milli>(high_resolution_clock::now() - t1).count() << " milliseconds \n";
			cout << "\nx:\n" << x;
			cout << "2-norm: " << norm(A * x - b)<<endl;
		}break;
		case 7: {
			int iters;
			cout << "number of iterations?: ";
			cin >> iters;
			t1 = high_resolution_clock::now();
			x = GaussSeidel(A, b, b, iters);
			cout << duration<double, std::milli>(high_resolution_clock::now() - t1).count() << " milliseconds \n";
			cout << "\nx:\n" << x;
			cout << "2-norm: " << norm(A * x - b) << endl;
		}break;
		case 8: {
			int iters;
			cout << "number of iterations?: ";
			cin >> iters;
			t1 = high_resolution_clock::now();
			x = SOR(A, b, b,1.2,iters);
			cout << duration<double, std::milli>(high_resolution_clock::now() - t1).count() << " milliseconds \n";
			cout << "\nx:\n" << x;
			cout << "2-norm: " << norm(A * x - b) << endl;
		}break;
		case 9: {
			t1 = high_resolution_clock::now();
			PowerIteration(A, eiganVal, eiganvec);
			cout << duration<double, std::milli>(high_resolution_clock::now() - t1).count() << " milliseconds \n";
			cout << "dominant eigan value: " << eiganVal<<endl;
			cout << "eigan vector :\n" << eiganvec << endl;
		}break;
		case 10: {
			int n;
			cout << "\nplease enter number of coefficeints: ";
			cin >> n;
			vector<term> arrf(n);
			cout << "for each term enter the coefficient then exponent \n(ex: \'1 2 -3 4\' -> \'1x^2-3x^4=0)\n";
			for (int i = 0; i < n; i++)
				cin >> arrf[i].c >> arrf[i].e;
			polynomial f(arrf);
			double guess;
			cout << "please enter a guess for the root: ";
			cin >> guess;
			double root = NewtonsMethod(f,guess);
			cout <<f<<endl;
			cout << "where x = " << root<<endl;
		}break;
		case 11: {
			int n;
			cout << "\nplease enter number of coefficeints: ";
			cin >> n;
			vector<term> arrf(n);
			cout << "for each term enter the coefficient then exponent \n(ex: \'1 2 -3 4\' -> \'1x^2-3x^4=0)\n";
			for (int i = 0; i < n; i++)
				cin >> arrf[i].c >> arrf[i].e;
			polynomial f(arrf);
			double xL, xR;
			cout << "what is xL?" << endl;
			cin >> xL;
			cout << "what is xR?" << endl;
			cin >> xR;
			double root = Bisection(f, xL,xR);
			cout << f << endl;
			cout << "where x = " << root << endl;
		}break;
		case 12: {
			int nx=sqrt(x.rows()), ny= sqrt(x.rows());
			ofstream out;
			out.open("x_graph.dat");
			for(int y=0;y<ny;y++)
				for (int xx = 0; xx < nx; xx++){
					out <<xx<<" "<<y<<" "<<x[0][xx + y * nx]<<"\n";
				}
			out.close();
			out.open("b_graph.dat");
			for (int y = 0; y<ny; y++)
				for (int xx = 0; xx < nx; xx++) {
					out << xx << " " << y << " " << b[0][xx + y * nx] << "\n";
				}
			out.close();

			out.open("function.txt");
			out << "set hidden3d\n"<<
				   "set dgrid3d 50, 50 qnorm 2\n"<<
				   "splot 'x_graph.dat' with lines,'b_graph.dat' with lines\n";
			out.close();

			system("gnuplot.exe -p function.txt ");
		}break;
		case 13: {
			cout << "\nA: \n" << A << "b:\n" << b << endl;
		}break;
		case 14: {
			cout << "\nA[" << A.rows() << "][" << A.cols() << "]\n";
			cout << "b[" << b.rows() << "]\n";
			if(!x.isEmpty())
				cout << "x[" << x.rows() << "]\n";
			else cout << "no x calculated\n";
			if (eiganVal != -0xff) {
				cout << "eigan vector[" << eiganvec.rows() << "]\n";
				cout << "eigan value = " << eiganVal<<"\n";
			}
			else cout << "no eigan value or vector calculated\n";
		}break;
		}


	}
}




int main() {
	interactiveMenu();
	cout << "done!\n";
	getchar();
}

