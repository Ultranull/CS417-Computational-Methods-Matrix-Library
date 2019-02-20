#pragma once

#include<math.h>
#include<iostream>
#include<vector>
#include<random>

using namespace std;

struct term {
	double e;
	double c;

	term(double c,double e):e(e),c(c){}
	term():term(0,0){}

	double eval(double x) {
		return c * pow(x, e);
	}

};

ostream& operator<<(ostream &out, term t) {
	if (t.c == 0.)
		return out;
	out << " ";
	if (t.c > 0.)
		out << "+";
	if (t.c != 1.)
		out << t.c;
	if(t.e != 0.)
		out << "x";
	if (t.e != 1.&&t.e != 0.)
		out << "^";
	
		if (t.e > 1.)
			out << t.e;
		else if (t.e < 0.)
			out << "(" << t.e << ")";
	
	return out;
}


struct polynomial {
	vector<term> function;
	polynomial(vector<term>f):function(f){}
	polynomial(int n) {
		function = vector<term>(n, term());
	}
	polynomial() {}

	double operator()(double x) {
		double sum = 0;
		for (int i = 0; i < function.size(); i++)
			sum += function[i].eval(x);
		return sum;
	}

	int size() {
		return function.size();
	}

	term& operator[](int i) {
		return function[i];
	}
};

ostream& operator<<(ostream &out, polynomial f) {
	for (int i = 0; i < f.size(); i++)
		out << f[i];
	return out;
}

polynomial ddx(polynomial f) {
	vector<term> fp;
	for (int i = 0; i < f.size(); i++)
		fp.push_back(term(f[i].c*f[i].e, f[i].e - 1));
	return polynomial(fp);
}

polynomial randomPolynom(default_random_engine &gen, double min, double max, int n = 4, bool rounded = false) {
	uniform_real_distribution<double> dist(min, max);
	polynomial out(n);
	for (int i = 0; i < n; i++) {
		if (rounded)
			out[i] = term(round(dist(gen)), round(dist(gen)));
		else
			out[i] = term(dist(gen), round(dist(gen)));
	}
	return out;
}

double NewtonsMethod(polynomial f, double guess = 1., int max = 100, double epsilon = -15) {
	double x0 = guess, xnew, error = abs(f(x0)), tol = pow(10, epsilon);
	int ct = 0;
	polynomial dfdx = ddx(f);
	while (error > tol&&ct < max) {
		xnew = x0 - (f(x0) / dfdx(x0));
		error = abs(f(xnew));
		x0 = xnew;
		ct++;
	}
	return xnew;
}