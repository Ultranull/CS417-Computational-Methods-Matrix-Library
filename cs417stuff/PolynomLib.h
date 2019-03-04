#pragma once

#include <math.h>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

struct term {
	double e;
	double c;

	term(double c, double e) :e(e), c(c) {}
	term() :term(0, 0) {}

	double eval(double x) {
		return c * pow(x, e);
	}

};

ostream& operator<<(ostream &out, term t) {
	if (t.c == 0.)
		return out;
	out << " ";
	if (t.c < 0.)
		out << "-";
	if (t.c != 1.)
		out << abs(t.c);
	if (t.e != 0.)
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
	polynomial(vector<term>f) :function(f) {}
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
	double highestOrder() {
		double max = 0;
		for (int i = 0; i < function.size(); i++)
			if (function[i].e > max)
				max = function[i].e;
		return max;
	}

	void sortFunction() {
		sort(function.begin(), function.end(), [](term a, term b) {
			return a.e > b.e;
		});
	}

	int size() {
		return function.size();
	}

	term& operator[](int i) {
		return function[i];
	}
};

ostream& operator<<(ostream &out, polynomial f) {
	for (int i = 0; i < f.size(); i++) {
		if (f[i].c == 0.)
			continue;
		out << " ";
		if (f[i].c > 0. && i > 0)
				out << "+ ";
		if (f[i].c < 0.) {
			out << "-";
			if (i > 0) out << " ";
		}
		if (f[i].c != 1. && f[i].c != -1.)
			out << abs(f[i].c);
		if (f[i].c != 1.&&f[i].e != 0)
			out << "*";
		if (f[i].e != 0.)
			out << "x";
		if (f[i].e != 1.&&f[i].e != 0.)
			out << "**";
		if (f[i].e > 1.)
			out << f[i].e;
		else if (f[i].e < 0.)
			out << "(" << f[i].e << ")";
	}
	//out << " = 0";
	return out;
}

polynomial ddx(polynomial f) {
	vector<term> fp;
	for (int i = 0; i < f.size(); i++)
		fp.push_back(term(f[i].c*f[i].e, f[i].e - 1));
	return polynomial(fp);
}

polynomial randomPolynom(default_random_engine &gen, double minc, double maxc, double mine, double maxe, int n = 4, bool rounded = false) {
	uniform_real_distribution<double> distc(minc, maxc);
	uniform_real_distribution<double> diste(mine, maxe);
	polynomial out(n);
	for (int i = 0; i < n; i++) {
		if (rounded)
			out[i] = term(round(distc(gen)), round(diste(gen)));
		else
			out[i] = term(distc(gen), round(diste(gen)));
	}
	return out;
}

double NewtonsMethod(polynomial f, double guess = 1., int max = 100, double epsilon = -15) {
	double x0 = guess, xnew = 0, error = abs(f(x0)), tol = pow(10, epsilon);
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