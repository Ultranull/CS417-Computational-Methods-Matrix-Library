#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <initializer_list>
#include <random>
#include <ctime>

using namespace std;

namespace dml {
	template <size_t cols, size_t rows, class TYPE>
	struct mat {
		vector<vector<TYPE>> data;
		mat(TYPE t) {
			data = vector<vector<TYPE>>(cols, 
				vector<TYPE>(rows,t));
		}
		mat() :mat(TYPE()) {}
		mat(const mat<cols, rows, TYPE> &m) {
			this->data = m.data;
		}
		mat(vector<vector<TYPE>> d){
			mat<rows, cols, TYPE> t;
			t.data = d;
			data = t.transpose().data;
		}
		TYPE operator()(int i=0) {
			return data[i][i];
		}
		vector<TYPE>& operator[](int i) {
			return data[i];
		}
		mat<1, rows, TYPE> col(unsigned int ind) { 
			mat<1, rows, TYPE> out;
			out.data[ind] = data[0];
			return out; 
		}
		mat<cols, 1, TYPE> row(unsigned int ind) {
			mat<cols, 1, TYPE> out;
			for (int i = 0; i < cols; i++)
				out.data[i][0] = data[i][ind];
			return out;
		}
		mat<rows, cols, TYPE> transpose() {
			mat<rows, cols, TYPE> out;
			for (int r = 0; r < rows; r++)
				for (int c = 0; c < cols; c++)
						out.data[r][c] = data[c][r];
			return out;
		}

		void insertROW(int ind,mat<cols, 1, TYPE> r){
			for (int i = 0; i < cols; i++)
				data[i][ind] = r.data[i][0];
		}
		void insertCOL(int ind,mat<1, rows, TYPE> r) {
				data[ind] = r.data[0];
		}
		int numcols() {
			return cols;
		}
		int numrows() {
			return rows;
		}
		static mat<cols, rows, TYPE> identity() {
			mat<cols, rows, TYPE> out;
			for (int r = 0; r < rows; r++)
				for (int c = 0; c < cols; c++)
					if (r == c)
						out.data[c][r] = 1;
			return out;
		}
	};
	template <size_t cols, size_t rows, class TYPE>
	ostream& operator<<(ostream &out, mat<cols, rows, TYPE> &m) {
		for (int r = 0; r < rows; r++) {
			out << "[";
			for (int c = 0; c < cols; c++)
				out <<setw(10)<< m.data[c][r] << " ";
			out << "]\n";
		}
		return out;
	}
	template <size_t cols, size_t rows, size_t cols2, class TYPE>
	mat<cols2, rows, TYPE> operator*(mat<cols, rows, TYPE> a, mat<cols2, cols, TYPE> b) {
		mat<cols2, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols2; j++)
				for(int k=0;k<cols;k++)
					out.data[j][i] += a.data[k][i]*b.data[j][k];
		return out;
	}
	template<size_t cols,size_t rows,class TYPE>
	mat<cols, rows, TYPE> operator*(mat<1, 1, TYPE> s, mat<cols, rows, TYPE> a) {
		mat<cols, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
					out.data[j][i] = s(0)*a.data[j][i];
		return out;
	}
	template<size_t cols, size_t rows, class TYPE>
	mat<cols, rows, TYPE> operator/(mat<cols, rows, TYPE> a, mat<1, 1, TYPE> s) {
		mat<cols, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				out.data[j][i] = a.data[j][i]/ s(0);
		return out;
	}
	template<size_t cols, size_t rows, class TYPE>
	mat<cols, rows, TYPE> operator*(TYPE s, mat<cols, rows, TYPE> a) {
		mat<cols, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				out.data[j][i] = s*a.data[j][i];
		return out;
	}
	template<size_t cols, size_t rows, class TYPE>
	mat<cols, rows, TYPE> operator/(mat<cols, rows, TYPE> a, TYPE s) {
		mat<cols, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				out.data[j][i] = a.data[j][i] / s;
		return out;
	}

	template<size_t cols, size_t rows, class TYPE>
	mat<cols, rows, TYPE> operator+(mat<cols, rows, TYPE> s, mat<cols, rows, TYPE> a) {
		mat<cols, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				out.data[j][i] = s.data[j][i] +a.data[j][i];
		return out;
	}
	template<size_t cols, size_t rows, class TYPE>
	mat<cols, rows, TYPE> operator-(mat<cols, rows, TYPE> s, mat<cols, rows, TYPE> a) {
		mat<cols, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				out.data[j][i] = s.data[j][i] - a.data[j][i];
		return out;
	}
	template<int cols, int rows, class TYPE>
	void swapROW(mat<cols, rows, TYPE> &x, int a, int b) {
		auto s = x.row(a);
		auto temp = x.row(b);
		x.insertROW(b, s);
		x.insertROW(a, temp);
	}
	template<int cols, int rows,class TYPE>
	void upperTriangular(mat<cols, rows, TYPE> &A, mat<1, rows, TYPE> &b) {
		for (int k = 0; k < A.numcols(); k++) {
			int ind = 0;
			TYPE max = 0;
			for (int i = 0; i < A.numrows(); i++)
				if (abs(A[k][i]) > max) {
					ind = i;
					max = abs(A[k][i]);
				}
			swapROW(A, ind, k);
			swapROW(b, ind, k);

			TYPE Ak = A(k);
			auto rowk = A.row(k);
			rowk = rowk / Ak;
			A.insertROW(k, rowk);

			auto bk = b.row(k);
			bk = bk / Ak;
			b.insertROW(k, bk);

			for (int r = k + 1; r < A.numrows(); r++) {
				TYPE val = A[k][r];
				for (int j = k; j < A.numcols(); j++) {
					A[j][r] = A[j][r] - (val*A[j][k]);
				}
				b[0][r] = b[0][r] - val * b[0][k];
			}
		}
	}

	template<int cols, int rows,class TYPE>
	mat<1, rows, TYPE> backSolve(mat<cols, rows, TYPE> A, mat<1, rows, TYPE> b) {
		mat<1, rows, TYPE> x(0.f);
		x[0][x.numrows() - 1] = b[0][b.numrows() - 1];
		for (int r = A.numrows() - 2; r >= 0; r--) {
			TYPE sum = 0;
			for (int c = r + 1; c < A.numcols(); c++) {
				sum += A[c][r] * x[0][c];
			}
			x[0][r] = b[0][r] - sum;
		}
		return x;
	}
	template<int cols, int rows, class TYPE>
	mat<1, rows, TYPE> forwordSolve(mat<cols, rows, TYPE> A, mat<1, rows, TYPE> b) {
		mat<1, rows, TYPE> x(0.f);
		x[0][0] = b[0][0]/A[0][0];
		for (int r = 1; r < A.numrows(); r++) {
			TYPE sum = 0;
			for (int c = 0; c < r; c++) {
				sum += A[c][r] * x[0][c];
			}
			x[0][r] = (b[0][r] - sum)/A[r][r];
		}
		return x;
	}
	template<int cols, int rows, class TYPE>
	mat<1, rows, TYPE> GaussianElimination(mat<cols, rows, TYPE> A, mat<1, rows, TYPE> b) {
		upperTriangular(A, b);
		return backSolve(A, b);
	}
	template<int n, class TYPE>
	void decomposeMatrix(mat<n, n, TYPE> A, mat<n, n, TYPE> &L, mat<n, n, TYPE> &U) {
		int i, j, k;
		double sum = 0;

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
	template<int cols, int rows, class TYPE>
	mat<cols, rows, TYPE> singularMatrix(TYPE min,TYPE max,int n = 1,bool rounded=false) {
		default_random_engine gen(time(NULL));
		uniform_real_distribution<TYPE> dist(min,max);
		mat<cols, rows, TYPE> out;

		for (int c = 0; c < cols; c++) {
			for (int r = 0; r < rows; r++) {
				if(rounded)
					out[c][r] = round(dist(gen));
				else
					out[c][r] = dist(gen);
			}
		}

		for (int c = 0; c < cols; c++) {
			TYPE sum = TYPE();
			for (int r = 0; r < rows; r++) {
				sum += abs(out[c][r]);
			}
			out[c][c] = sum * n;
		}
		return out;
	}
	template<int rows, class TYPE>
	mat<1, rows, TYPE> randomVector(TYPE min, TYPE max, bool rounded = false) {
		default_random_engine gen(time(NULL));
		uniform_real_distribution<TYPE> dist(min, max);
		mat<1, rows, TYPE> out;
		for (int i = 0; i < rows; i++) {
			if (rounded)
				out[0][i] = round(dist(gen));
			else
				out[0][i] = (dist(gen));
		}
		return out;
	}
}