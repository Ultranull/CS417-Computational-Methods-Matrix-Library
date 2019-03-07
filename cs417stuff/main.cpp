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


#include "MatrixLib.h"
#include "PolynomLib.h"

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
	cout << "8. Solve current matrix with SOR method\n";
	cout << "9. Find EiganVector and dominant Eigan Value using power iteration\n";
	cout << "10. Find a root to polynomial using newtons method\n";
	//cout << "9. Save current solution vector x to file\n";
	//cout << "10. Graph current solution vector x using GNUPLOT\n";
	cout << "13. display A and b\n";
	cout << "0. exit\n";

	cin >> option;
	return option;
}
void interactiveMenu() {
	default_random_engine gen(time(NULL)); 
	high_resolution_clock::time_point t1;
	mat A(1,2,2), b(1, 1, 2),L,U,D,x,eiganvec;
	double eiganVal;
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
		case 13: {
			cout << "\nA: \n" << A << "b:\n" << b << endl;
		}break;
		}


	}
}




int main() {
	const size_t N = 2;
	double min = -5, max = 10;
	default_random_engine gen(time(NULL));
	{
		//const size_t N = 2;
		//double min = 0, max = 10;
		//default_random_engine gen(time(NULL));
		//if(0){
		//	
		//	mat A = nonsingularMatrix(gen, min, max, N, true);
		//	mat b = randomVector(gen, min, max, N, true);


		//	mat bk;
		//	double mu;

		//	PowerIteration(A, mu, bk);

		//	mat Al = A * bk, lv = mu * bk;
		//	cout << norm(Al) << endl << norm(lv) << endl;

		//	mat f = randomVector(gen, min, max, 100, true).transpose();

		//	double x = NewtonsMethod(f);
		//	cout << x << endl;
		//	cout << abs(fx(f, x)) << endl;
		//	/*for (int i = 0; i < 100; i++) {
		//		mat At = nonsingularMatrix(gen, min, max, N, true);
		//		mat bt = randomVector(gen, min, max, N, true);
		//		vector<double> error;
		//		mat xnew = SOR(At, bt, bt,.2, error, 1, 100);
		//		mat twoNorm = At * xnew - bt;
		//		cout << "SOR: " <<error.size()<<" 2norm: "<<norm(twoNorm) << endl;
		//		error.clear();
		//		xnew = GaussSeidel(At, bt, bt, error, 1, 100);
		//		twoNorm = At * xnew - bt;
		//		cout <<"guass-siedel: "<< error.size() << " 2norm: " << norm(twoNorm) << endl << endl;
		//	}*/
		//}
	}
	//interactiveMenu();
	cout << endl;
	{

		polynomial f = randomPolynom(gen, -5, 5,0,3, 8,true);
		normal_distribution<double> nd(0, 1);
		normal_distribution<double> nd0(0., 15.);
		uniform_real_distribution<double> dist(0, 3);
		vector<vector<double>> Amat;
		vector<vector<double>> bvec;

		ofstream out;
		out.open("data.txt");
		double num = 5;
		for (double x = -num/2; x < num / 2; x+=.01) {
			for (int i = 0; i < 1; i++) {
				double rand0 = nd0(gen);
				double rand =  nd(gen);
				double xp = rand + x;
				double b = rand0 + f(xp);
				cout << rand0 << endl;
				out << xp << " " << b << endl;
				vector<double> Arow;
				for (int i = 0; i < f.highestOrder() + 1; i++)
					Arow.push_back(pow(xp, i));
				Amat.push_back(Arow);
				bvec.push_back(vector<double>(1, b));
			}
		}
		out.close();
		mat A(Amat);
		mat b(bvec);

		mat At = A.transpose();
		mat AtA = At * A;
		mat Atb = At * b;

		mat x = GaussianElimination(AtA, Atb);

		vector<term> fpv;
		for (int i = x.rows()-1; i >=0; i--)
			fpv.push_back(term(x[0][i], i));
		polynomial fp(fpv);
		
		f.sortFunction();
		cout << f.highestOrder() << endl;
		cout << f << endl;
		cout << fp << endl;
		


		out.open("function.txt");
		out << "fp(x)="<< fp<<endl;
		out << "f(x)=" << f << endl;
		out << "plot 'data.txt' with points,\\\n";
		out << "     fp(x) lt rgb \"red\",\\\n";
		out << "     f(x)\n";
		out.close();

		//system("gnuplot.exe -p function.txt ");
		
	}


	cout << "done!\n";
	getchar();
}

