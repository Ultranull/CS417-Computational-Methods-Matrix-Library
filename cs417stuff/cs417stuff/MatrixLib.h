#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <initializer_list>

using namespace std;

namespace dml {
/*
	template <size_t cols, class TYPE>
	struct vec {
		vector<TYPE> data;
		vec() :vec(TYPE()) {
		}
		vec(TYPE t) {
			data = vector<TYPE>(cols, t);
		}
		vec(initializer_list<TYPE> l) {
			for (int i = 0; i < cols; i++)
				data.push_back(*(l.begin() + i));
		}
		vec(const vec<cols, TYPE> &a) {
			this->data = a.data;
		}
		TYPE& operator[](unsigned int i) { return data[i]; }
	};
	template <size_t cols, class TYPE>
	vec<cols, TYPE> operator+(vec<cols, TYPE> a, vec<cols, TYPE> b) {
		vec<cols, TYPE> out;
		for (int i = 0; i < cols; i++)
			out[i] = a[i] + b[i];
		return out;
	}
	template <size_t cols, class TYPE>
	vec<cols, TYPE> operator-(vec<cols, TYPE> a, vec<cols, TYPE> b) {
		vec<cols, TYPE> out;
		for (int i = 0; i < cols; i++)
			out[i] = a[i] - b[i];
		return out;
	}
	template <size_t cols, class TYPE>
	vec<cols, TYPE> operator*(TYPE a, vec<cols, TYPE> b) {
		vec<cols, TYPE> out;
		for (int i = 0; i < cols; i++)
			out[i] = a * b[i];
		return out;
	}
	template <size_t cols, class TYPE>
	TYPE dot(vec<cols, TYPE> a, vec<cols, TYPE> b) {
		TYPE out = TYPE();
		for (int i = 0; i < cols; i++)
			out += a[i] * b[i];
		return out;
	}
	template <size_t cols, class TYPE>
	ostream& operator<<(ostream &out, vec<cols, TYPE> &v) {
		for (int i = 0; i < cols; i++)
			out << v[i] << " ";
		return out;
	}*/
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
		TYPE operator()(void) {
			return data[0][0];
		}
		mat<1, rows, TYPE> col(unsigned int ind) { 
			mat<1, rows, TYPE> out;
			for (int i = 0; i < rows; i++)
				out.data[ind][i] = data[0][i];
			return out; 
		}
		mat<cols, 1, TYPE> row(unsigned int ind) {
			mat<cols, 1, TYPE> out;
			for (int i = 0; i < cols; i++)
				out.data[i][ind] = data[i][0];
			return out;
		}
		mat<rows, cols, TYPE> transpose() {
			mat<rows, cols, TYPE> out;
			for (int r = 0; r < rows; r++)
				for (int c = 0; c < cols; c++)
						out.data[r][c] = data[c][r];
			return out;
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
				out <<setw(3)<< m.data[c][r] << " ";
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
	mat<cols, rows, TYPE> operator*(TYPE s, mat<cols, rows, TYPE> a) {
		mat<cols, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
					out.data[j][i] = s*a.data[j][i];
		return out;
	}

}