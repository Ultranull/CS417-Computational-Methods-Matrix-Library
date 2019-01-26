#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <initializer_list>

#include "MatrixLib.h"

using namespace std;

namespace veclib {
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
	}
	template <size_t cols, size_t rows, class TYPE>
	struct mat {
		vector<vec<cols, TYPE>> data;
		mat(TYPE t) {
			data = vector<vec<cols, TYPE>>(cols, vec<cols, TYPE>(t));
		}
		mat() :mat(TYPE()) {}
		mat(initializer_list<vec<cols, TYPE>> l) {
			for (int i = 0; i < rows; i++)
				data.push_back(*(l.begin() + i));
		}
		mat(const mat<cols, rows, TYPE> &m) {
			this->data = m.data;
		}
		vec<cols, TYPE>& operator[](unsigned int i) { return data[i]; }
		vec<rows, TYPE> col(unsigned int ind) {
			vec<rows, TYPE> out;
			for (int i = 0; i < rows; i++)
				out[i] = data[i][ind];
			return out;
		}
		static mat<cols, rows, TYPE> identity() {
			mat<cols, rows, TYPE> out;
			for (int r = 0; r < rows; r++)
				for (int c = 0; c < cols; c++)
					if (r == c)
						out[r][c] = 1;
			return out;
		}
	};
	template <size_t cols, size_t rows, class TYPE>
	ostream& operator<<(ostream &out, mat<cols, rows, TYPE> &m) {
		for (int i = 0; i < rows; i++)
			out << m[i] << "\n";
		return out;
	}
	template <size_t cols, size_t rows, class TYPE>
	vec<rows, TYPE> operator*(mat<cols, rows, TYPE> a, vec<cols, TYPE> b) {
		vec<rows, TYPE> out;
		for (int i = 0; i < cols; i++)
			out[i] = dot(a[i], b);
		return out;
	}
	template <size_t cols, size_t rows, class TYPE>
	mat<cols, rows, TYPE> operator*(mat<cols, rows, TYPE> a, mat<cols, rows, TYPE> b) {
		mat<cols, rows, TYPE> out;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				out[i][j] = dot(a[i], b.col(j));
		return out;
	}
	typedef mat<4, 4, float> mat4;
	typedef vec<4, float> vec4;

}
int main() {
	{
		veclib::vec4 z{ 1,2,3,1 };
		veclib::mat4 x{
			{1,0,0,0},
			{0,1,0,5},
			{0,0,1,0},
			{0,0,0,1},
		};
		veclib::mat4 y{
			{3,0,0,0},
			{0,1,0,0},
			{0,0,2,0},
			{0,0,0,1},
		};
		auto out = x * y * z;
		cout << out << endl;
		cout << endl;
	} {
		using namespace dml;
		typedef mat<4, 4, float> mat4;
		typedef mat<1, 4, float> scalor;

		mat4 x = mat4::identity();
		scalor y(2);
		auto z = x * y;
		cout << z;

	}
	//this is a test to see if i can get around not using git
	getchar();
}

