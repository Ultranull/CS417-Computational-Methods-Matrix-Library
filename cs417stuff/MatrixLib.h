#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <ctime>

using namespace std;

namespace dml {
	struct mat {
		vector<vector<double>> data;
		mat(double t,int cols,int rows) {
			data = vector<vector<double>>(cols,
				vector<double>(rows,t));
		}
		mat(){}
		mat(int cols,int rows):mat(0,cols,rows){}
		mat(const mat &m) {
			this->data = m.data;
		}
		mat(vector<vector<double>> d){
			mat t(0,d.size(),d[0].size());
			t.data = d;
			data = t.transpose().data;
		}
		double operator()(int i=0) {
			return data[i][i];
		}
		vector<double>& operator[](int i) {
			return data[i];
		}
		mat col(unsigned int ind) { 
			mat out(0,1,rows());
			out.data[ind] = data[0];
			return out; 
		}
		mat row(unsigned int ind) {
			mat out(0,cols(),1);
			for (int i = 0; i < cols(); i++)
				out.data[i][0] = data[i][ind];
			return out;
		}
		mat transpose() {
			mat out(0,rows(), cols());
			for (int r = 0; r < rows(); r++)
				for (int c = 0; c < cols(); c++)
						out.data[r][c] = data[c][r];
			return out;
		}

		void insertROW(int ind,mat r){
			for (int i = 0; i < cols(); i++)
				data[i][ind] = r.data[i][0];
		}
		void insertCOL(int ind,mat r) {
				data[ind] = r.data[0];
		}
		int cols() {
			return data.size();
		}
		int rows() {
			return data[0].size();
		}
		static mat identity(int n) {
			mat out(0,n,n);
			for (int r = 0; r < n; r++)
				for (int c = 0; c < n; c++)
					if (r == c)
						out.data[c][r] = 1;
			return out;
		}
	};
	ostream& operator<<(ostream &out, mat &m) {
		for (int r = 0; r < m.rows(); r++) {
			out << "[";
			for (int c = 0; c < m.cols(); c++)
				out <<setw(10)<< m.data[c][r] << " ";
			out << "]\n";
		}
		return out;
	}

	double norm(mat v) {
		double sum = 0;
		for (int i = 0; i < v.rows(); i++) {
			sum += v[0][i] * v[0][i];
		}
		return sqrt(sum);
	}

	mat operator*(mat a, mat b) {
		mat out(0,b.cols(),a.rows());
		for (int i = 0; i < a.rows(); i++)
			for (int j = 0; j < b.cols(); j++)
				for(int k=0;k<a.cols();k++)
					out.data[j][i] += a.data[k][i]*b.data[j][k];
		return out;
	}

	mat operator*(double s, mat a) {
		mat out = a;
		for (int i = 0; i < a.rows(); i++)
			for (int j = 0; j < a.cols(); j++)
				out.data[j][i] = s*a.data[j][i];
		return out;
	}
	mat operator/(mat a, double s) {
		mat out=a;
		for (int i = 0; i < a.rows(); i++)
			for (int j = 0; j < a.cols(); j++)
				out.data[j][i] = a.data[j][i] / s;
		return out;
	}

	mat operator/(double s,mat a) {
		mat out=a;
		for (int i = 0; i < a.rows(); i++)
			for (int j = 0; j < a.cols(); j++) {
				if(a.data[j][i]==0)
					out.data[j][i] = 0;
				else
					out.data[j][i] = s / a.data[j][i];
			}
		return out;
	}

	mat operator+(mat s, mat a) {
		mat out=a;
		for (int i = 0; i < a.rows(); i++)
			for (int j = 0; j < a.cols(); j++)
				out.data[j][i] = s.data[j][i] +a.data[j][i];
		return out;
	}
	mat operator-(mat s, mat a) {
		mat out=s;
		for (int i = 0; i < a.rows(); i++)
			for (int j = 0; j < a.cols(); j++)
				out.data[j][i] = s.data[j][i] - a.data[j][i];
		return out;
	}
	void swapROW(mat &x, int a, int b) {
		mat s = x.row(a);
		mat temp = x.row(b);
		x.insertROW(b, s);
		x.insertROW(a, temp);
	}
	void upperTriangular(mat &A, mat &b) {
		for (int k = 0; k < A.cols(); k++) {
			int ind = 0;
			double max = 0;
			for (int i = 0; i < A.rows(); i++)
				if (abs(A[k][i]) > max) {
					ind = i;
					max = abs(A[k][i]);
				}
			swapROW(A, ind, k);
			swapROW(b, ind, k);


			double Ak = A(k);
			mat rowk = A.row(k);
			rowk = rowk / Ak;
			A.insertROW(k, rowk);

			mat bk = b.row(k);
			bk = bk / Ak;
			b.insertROW(k, bk);

			for (int r = k + 1; r < A.rows(); r++) {
				double val = A[k][r];
				for (int j = k; j < A.cols(); j++) {
					A[j][r] = A[j][r] - (val*A[j][k]);
				}
				b[0][r] = b[0][r] - val * b[0][k];
			}
		}
	}

	mat backSolve(mat A, mat b) {
		mat x=b;
		x[0][x.rows() - 1] = b[0][b.rows() - 1];
		for (int r = A.rows() - 2; r >= 0; r--) {
			double sum = 0;
			for (int c = r + 1; c < A.cols(); c++) {
				sum += A[c][r] * x[0][c];
			}
			x[0][r] = b[0][r] - sum;
		}
		return x;
	}

	mat forwordSolve(mat A, mat b) {
		mat x=b;
		x[0][0] = b[0][0]/A[0][0];
		for (int r = 1; r < A.rows(); r++) {
			double sum = 0;
			for (int c = 0; c < r; c++) {
				sum += A[c][r] * x[0][c];
			}
			x[0][r] = (b[0][r] - sum)/A[r][r];
		}
		return x;
	}

	mat GaussianElimination(mat A, mat b) {
		upperTriangular(A, b);
		return backSolve(A, b);
	}

	void decomposeMatrix(mat A, mat &L, mat &U) {
		int i, j, k;
		double sum = 0;
		int n=A.rows();
		L = mat(n, n);
		U = mat(n, n);
		for (i = 0; i < n; i++) {
			U[i][i] = 1;
		}

		for (j = 0; j < n; j++) {
			for (i = j; i < n; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L[k][i] * U[j][k];
				}
				L[j][i] = A[j][i] - sum;
			}

			for (i = j; i < n; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L[k][j] * U[i][k];
				}
				U[i][j] = (A[i][j] - sum) / L[j][j];
			}
		}
	}
	mat LUDecomposition(mat A, mat b) {
		mat L, U;
		decomposeMatrix(A, L, U);
		mat y = forwordSolve(L, b);;
		return backSolve(U, y);
	}
	mat nonsingularMatrix(default_random_engine &gen,double min,double max,int n = 4,bool rounded=false) {
		uniform_real_distribution<double> dist(min,max);
		mat out(0,n,n);
		for (int c = 0; c < n; c++) {
			for (int r = 0; r < n; r++) {
				if(rounded)
					out[c][r] = round(dist(gen));
				else
					out[c][r] = dist(gen);
			}
		}

		for (int c = 0; c < n; c++) {
			double sum = double();
			for (int r = 0; r < n; r++) {
				sum += abs(out[c][r]);
			}
			out[c][c] = sum;
		}
		return out;
	}

	mat randomVector(default_random_engine &gen,double min, double max, int n = 4, bool rounded = false) {
		uniform_real_distribution<double> dist(min, max);
		mat out(0,1,n);
		for (int i = 0; i < n; i++) {
			if (rounded)
				out[0][i] = round(dist(gen));
			else
				out[0][i] = (dist(gen));
		}
		return out;
	}

	void split(mat A, mat &L, mat &D, mat &U) {
		int n = A.cols();
		L = mat(n, n);
		U = mat(n, n);
		D = mat(n, n);
		for (int r = 0; r < A.rows(); r++) {
			for (int c = 0; c < A.cols(); c++) {
				if (r == c)
					D[c][r] = A[c][r];
				else if (r < c)
					U[c][r] = A[c][r];
				else if (r > c)
					L[c][r] = A[c][r];
			}
		}
	}

	mat JacobiIterative(mat A, mat b, mat guess, int iters = 100) {
		int N = A.cols();
		mat L(N, N), D(N, N), U(N, N);
		split(A, L, D, U);
		mat Dinv = 1. / D, LU = L + U;
		mat xold = guess, xnew(1, A.rows());
		unsigned int i = 0;
		while (abs(norm(xnew) - norm(xold)) > pow(10, -50) && i < iters) {
			xold = xnew;
			xnew = Dinv * (b - LU * xold);
			i++;
		}
		return xold;
	}

	mat JacobiIterative(mat A, mat b, mat guess,vector<double> &twoN,int errordiv,int iters = 100) {
		mat L, D, U;
		split(A, L, D, U);
		mat Dinv = 1. / D, LU = L + U;
		mat xold = guess, xnew(1,A.rows());
		unsigned int i = 0;
		while (abs(norm(xnew) - norm(xold)) > pow(10, -50) && i<iters) {
			xold = xnew;
			xnew = Dinv * (b - LU * xold);
			if (i%errordiv == 0)
				twoN.push_back(norm(A * xnew - b));
			i++;
		} 
		return xnew;
	}


	mat GaussSeidel(mat A, mat b, mat guess, vector<double> &twoN, int errordiv, int iters = 100) {
		int n = A.cols();
		mat xold = guess, xnew(0, 1, n);
		unsigned int c = 0;
		while (abs(norm(xnew) - norm(xold)) > pow(10, -50) && c < iters) {
			xnew = xold;
			for (int i = 0; i < n; i++) {
				double theta=0;
				for (int j = 0; j <n; j++){
					if (j != i)
						theta = theta + A[j][i] * xold[0][j];
				}
				xold[0][i] = (1 / A[i][i])*(b[0][i] - theta);
			}
			if (c%errordiv == 0)
				twoN.push_back(norm(A * xold - b));
			c++;
		}
		return xold;
	}

	mat GaussSeidel(mat A, mat b, mat guess, int iters = 100) {
		int n = A.cols();
		mat xold = guess, xnew(0, 1, n);
		unsigned int c = 0;
		while (abs(norm(xnew) - norm(xold)) > pow(10, -50) && c < iters) {
			xnew = xold;
			for (int i = 0; i < n; i++) {
				double theta = 0;
				for (int j = 0; j < n; j++) {
					if (j != i)
						theta = theta + A[j][i] * xold[0][j];
				}
				xold[0][i] = (1 / A[i][i])*(b[0][i] - theta);
			}
			c++;
		}
		return xold;
	}
	mat SOR(mat A, mat b, mat guess,double con, vector<double> &twoN, int errordiv, int iters = 100) {
		int n = A.cols();
		mat xold = guess, xnew(0, 1, n);
		unsigned int c = 0;
		double omega=0;
		while (abs(norm(xnew) - norm(xold)) > pow(10, -50)&&c < iters) {
			xnew = xold;
			for (int i = 0; i < n; i++) {
				double theta = 0;
				for (int j = 0; j < n; j++) { 
					if (j != i)
						theta = theta + A[j][i] * xold[0][j];
				}
				omega = con;
				xold[0][i] = (1. - omega)*xold[0][i]+(omega / A[i][i])*(b[0][i] - theta);
			}
			if (c%errordiv == 0)
				twoN.push_back(norm(A * xold - b));
			c++;
		}
		return xold;
	}
}