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
#include <sstream>
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
	cout << "13. Save current solution vector x using GNUPLOT\n";
	cout << "14. display A and b\n";
	cout << "15. display stored data information\n";
	cout << "16. create log file\n";
	cout << "17. load a command file\n";
	cout << "0. exit\n";

	cin >> option;
	return option;
}
void interactiveMenu() {
	ofstream log;
	bool makelog = false;
	default_random_engine gen(time(NULL)); 
	high_resolution_clock::time_point t1;
	mat A(1,2,2), b(1, 1, 2),L,U,D,x,eiganvec;
	double eiganVal=-0xff;
	string fn="out";
	ifstream commandfile; 
	streambuf *backup=cin.rdbuf();
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
			cout << "what file?: ";
			cin >> fn;
			load(A, b, fn);
			if (makelog)log << "\nopening file:" << fn << endl;
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
			if (makelog)log << "saving file:" << n << endl;

		}break;
		case 4: {
			t1 = high_resolution_clock::now();
			x = GaussianElimination(A, b);
			double time = duration<double, std::milli>(high_resolution_clock::now() - t1).count();
			cout << time << " milliseconds \n";
			double norm2 = norm(A * x - b);
			cout << "2-norm: " << norm2 << endl;
			if (makelog)log << "GaussianElimination in " << time << " milliseconds \n";
		}break;
		case 5: {
			t1 = high_resolution_clock::now();
			x = LUDecomposition(A, b);
			double time = duration<double, std::milli>(high_resolution_clock::now() - t1).count();
			cout << time << " milliseconds \n";
			if (makelog)log << "LUDecomposition in " << time << " milliseconds \n";
		}break;
		case 6: {
			int iters,count;
			cout << "number of iterations?: ";
			cin >> iters;
			t1 = high_resolution_clock::now();
			x = JacobiIterative(A, b, b,count , iters);
			double time = duration<double, std::milli>(high_resolution_clock::now() - t1).count();
			cout << time << " milliseconds \n";
			double norm2 = norm(A * x - b);
			cout << "2-norm: " << norm2<<endl;
			cout << "number of iterations: " << count << "\n";
			if (makelog)log << "JacobiIterative in " << time << " milliseconds \n" <<
							   "2norm = " << norm2 << endl<<
							   "iterations = " << count << endl;
		}break;
		case 7: {
			int iters,count;
			cout << "number of iterations?: ";
			cin >> iters;
			t1 = high_resolution_clock::now();
			x = GaussSeidel(A, b, b,count , iters);
			double time = duration<double, std::milli>(high_resolution_clock::now() - t1).count();
			cout << time << " milliseconds \n";
			double norm2 = norm(A * x - b);
			cout << "2-norm: " << norm2 << endl;
			cout << "number of iterations: " << count << "\n";
			if (makelog)log << "GaussSeidel in " << time << " milliseconds \n" <<
				"2norm = " << norm2 << endl <<
				"iterations = " << count << endl;
		}break;
		case 8: {
			int iters;
			cout << "number of iterations?: ";
			cin >> iters;
			double omega;
			cout << "omega?: ";
			cin >> omega;
			int count;
			t1 = high_resolution_clock::now();
			x = SOR(A, b, b, omega,count,iters);
			double time = duration<double, std::milli>(high_resolution_clock::now() - t1).count();
			cout << time << " milliseconds \n";
			double norm2 = norm(A * x - b);
			cout << "2-norm: " << norm2 << endl;
			cout << "number of iterations: " << count << "\n";
			if (makelog)log << "SOR in " << time << " milliseconds \n" <<
				"2norm = " << norm2 << endl <<
				"omega = " << omega << endl <<
				"iterations = " << count << endl;
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
			double root = 0;
			if (makelog) {
				log << "\nnewtons method with " << f <<" = 0"<< endl;
				root = NewtonsMethod(log,f, guess,100,-5);
			}else
				root = NewtonsMethod(f, guess);
			cout <<f<<endl;
			cout << "where x = " << root<<endl;
			if (makelog)log << "where x = " << root << endl;
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
			int nx = sqrt(x.rows()), ny = sqrt(x.rows());
			ofstream out, outexcel;
			out.open(fn + "_x_gnuplot.txt");
			outexcel.open(fn + "_x_excel.txt");
			float min = 0xff, max = -0xff;
			for (int xx = 0; xx < nx; xx++) {
				for (int y = 0; y < ny; y++) {
					float z = x[0][xx + y * nx];
					out << xx << " " << y << " " << z << "\n";
					outexcel << z << " ";
					min = z < min ? z : min;
					max = z > max ? z : max;
				}
				outexcel << "\n";
				out << "\n";
			}
			out.close();
			outexcel.close();
			out.open("function.txt");
			out << //"set hidden3d\n"<<
				"set pm3d\n" <<
				"set palette defined (" << min << " \"blue\", " << (min + max) / 2 << " \"green\", " << max << " \"red\")\n" <<
				"set pal maxcolors 50\n" <<
				"splot '" + fn + "_x_gnuplot.txt' with pm3d\n";
			out.close();

			system("gnuplot.exe -p function.txt ");
		}break;
		case 13: {
			int nx = sqrt(x.rows()), ny = sqrt(x.rows());
			ofstream out,outexcel;
			out.open(fn+"_x_gnuplot.txt");
			outexcel.open(fn + "_x_excel.txt");
			float min = 0xff, max = -0xff;
			for (int xx = 0; xx < nx; xx++) {
				for (int y = 0; y < ny; y++) {
					float z = x[0][xx + y * nx];
					out << xx << " " << y << " " << z << "\n";
					outexcel << z<<" ";
					min = z < min ? z : min;
					max = z > max ? z : max;
				}
				outexcel << "\n";
				out << "\n";
			}
			out.close(); 
			outexcel.close();
			if (makelog)log << "saving image to " << fn << ".png\n";
			out.open("function.txt");
			out << //"set hidden3d\n"<<
				"set term png\n"<<
				"set output \""<<fn<<".png\"\n"<<
				"set pm3d\n" <<
				"set palette defined (" << min << " \"blue\", " << (min + max) / 2 << " \"green\", " << max << " \"red\")\n" <<
				"set pal maxcolors 50\n" <<
				"splot '"+ fn + "_x_gnuplot.txt' with pm3d\n";
			out.close();

			system("gnuplot.exe -p function.txt ");
		}break;
		case 14: {
			cout << "\nA: \n" << A << "b:\n" << b << endl;
			if (makelog)log << "\nA: \n" << A << "b:\n" << b << endl;
		}break;
		case 15: {
			cout << "\nA[" << A.rows() << "][" << A.cols() << "]\n";
			cout << "b[" << b.rows() << "]\n";
			if(!x.isEmpty())
				cout << "x[" << x.rows() << "]\n";
			else cout << "no x calculated\n";
			if (makelog) {
				log << "A[" << A.rows() << "][" << A.cols() << "]\n";
				log << "b[" << b.rows() << "]\n";
				if (!x.isEmpty())
					log << "x[" << x.rows() << "]\n";
				else log << "no x calculated\n";
			}
			if (eiganVal != -0xff) {
				cout << "eigan vector[" << eiganvec.rows() << "]\n";
				cout << "eigan value = " << eiganVal<<"\n";
			}
			else cout << "no eigan value or vector calculated\n";
		}break;
		case 16: {
			stringstream ss;
			double time = duration<double, std::nano>(high_resolution_clock::now() - t1).count();
			ss << "log" << (time) << ".txt";
			log.open(ss.str());
			makelog = true;
		}break;
		case 17: {
			string n;
			cout << "command file name: \n";
			cin >> n;
			commandfile.open(n);
			if (!commandfile) {
				cout << "error: command file not found!\n";
			}
			backup = cin.rdbuf();
			cin.rdbuf(commandfile.rdbuf());
		}break;
		case 100: {
			if (makelog) {
				log << "\nbistection of x + ln(x)=0\n";
				double xx=Bisection([](double x)->double {return std::log(x) + x; }, .1, 1.,100,log);
				log << "solution x = " << xx << endl;
			}
		}break;
		}
		if (cin.peek() == -1) {
			cin.rdbuf(backup);
			cout << "command file completed\a\n";
		}
	}
	commandfile.close();
	log.close();
}




int main() {
	interactiveMenu();

	/*ifstream fin;
	fin.open("x_data.txt", ios::in);
	int n;
	fin >> n;
	mat x = mat(1, n);
	double val;
	for (int i = 0; i < n; i++) {
		fin >> val;
		x[0][i] = val;
	}
	fin.close();

	default_random_engine gen(time(NULL));
	mat A = nonsingularMatrix(gen, 0, 1, n);
	mat b = A * x;

	save(A, b, "test2.txt");*/

	cout << "done!\n";
	//getchar();
}

