#define _CRT_SECURE_NO_WARNINGS
#include <cstdlib>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include "Matrix.h"

using namespace std;

vector< vector<double> > M;

const double eps = 1.e-10;


int main(int argc, char** argv) {
	Matrix m(5);
	for (int i = 0; i < 5; ++i){
		for (int j = 0; j < 5; ++j){
			double c = 0.1 * 23 * (i + 1) * (j + 1);
			c += 1;
			m.m[i][j] = 11.7 / (c * c * c * c * c * c * c);
		}
	}

	for (int i = 0; i < 5; ++i){
		for (int j = 0; j < 5; ++j)
			cout << scientific << m.m[i][j] << " ";
		cout << endl;
	}
	cout << endl;
	vector<double> f(5, 5);
	pair<vector<double>, Matrix> r = m.Gauss(f);
	cout << endl << "Condition number\n" << r.second.Norm()*m.Norm() << endl << endl << "A^(-1)\n";
	for (int i = 0; i < 5; ++i){
		for (int j = 0; j < 5; ++j)
			cout << scientific << r.second.m[i][j] << " ";
		cout << endl;
	}


	cout << "res" << endl;
	for (int i = 0; i < 5; ++i){
		cout << r.first[i] << " ";
	}
	cout << endl;
	cout << "det\n" << m.det() << endl;
	cout << endl << "A*A^(-1)\n";
	Matrix t(r.second * m);
	for (int i = 0; i < 5; ++i){
		for (int j = 0; j < 5; ++j)
			cout << t.m[i][j] << " ";
		cout << endl;
	}
	cout << endl;
	cout << endl << "Nav\n";
	vector<double> d(m * r.first);
	for (int i = 0; i < 5; ++i)
		cout << d[i] - f[i] << " ";
	cout << endl;
	cout << endl << endl;
	Matrix E(5);
	E.SetIdentity();
	Matrix tmp(m + E);
	Matrix zb(m + E * 1e-25);
	vector<double> z(zb.Gauss(f).first);
	cout << "zb" << endl;
	for (int i = 0; i < 5; ++i)
		cout << r.first[i] - z[i] << " ";
	cout << endl;
	vector<double> j(tmp.Jacobi(f));
	cout << "A + E" << endl;
	for (int i = 0; i < 5; ++i){
		for (int j = 0; j < 5; ++j)
			cout << scientific << tmp.m[i][j] << " ";
		cout << endl;
	}
	cout << "Jacobi res\n";
	for (int i = 0; i < 5; ++i){
		cout << j[i] << ' ';
	}
	cout << endl << "Nav\n";
	for (int i = 0; i < 5; ++i)
		cout << j[i] - f[i] << " ";
	cout << endl;
	return 0;
}

