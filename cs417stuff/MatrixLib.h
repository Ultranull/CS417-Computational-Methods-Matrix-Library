#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <initializer_list>

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
	template<int cols, int rows>
	void upperTriangular(mat<cols, rows, float> &A, mat<1, rows, float> &b) {
		for (int k = 0; k < A.numcols(); k++) {
			int ind = 0;
			float max = 0;
			for (int i = 0; i < A.numrows(); i++)
				if (abs(A[k][i]) > max) {
					ind = i;
					max = abs(A[k][i]);
				}
			swapROW(A, ind, k);
			swapROW(b, ind, k);

			float Ak = A(k);
			auto rowk = A.row(k);
			rowk = rowk / Ak;
			A.insertROW(k, rowk);

			auto bk = b.row(k);
			bk = bk / Ak;
			b.insertROW(k, bk);

			for (int r = k + 1; r < A.numrows(); r++) {
				float val = A[k][r];
				for (int j = k; j < A.numcols(); j++) {
					A.data[j][r] = A[j][r] - (val*A[j][k]);
				}
				b.data[0][r] = b[0][r] - val * b[0][k];
			}
		}
	}

	template<int cols, int rows>
	mat<1, rows, float> backSolve(mat<cols, rows, float> A, mat<1, rows, float> b) {
		mat<1, rows, float> x(0.f);
		x.data[0][x.numrows() - 1] = b[0][b.numrows() - 1];
		for (int r = A.numrows() - 2; r >= 0; r--) {
			float sum = 0;
			for (int c = r + 1; c < A.numcols(); c++) {
				sum += A[c][r] * x[0][c];
			}
			x.data[0][r] = b[0][r] - sum;
		}
		return x;
	}
}