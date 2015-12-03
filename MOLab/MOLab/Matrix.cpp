#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>

const double eps = 1.e-30;

Matrix::Matrix() {
	Det = nan("");
}

Matrix::Matrix(int n) {
	Det = 0;
	m.resize(n, std::vector<double>(n, 0));
}

Matrix::Matrix(const Matrix& orig) {
	Det = orig.Det;
	m.assign(orig.m.begin(), orig.m.end());
}

void Matrix::SetIdentity() {
	SetNull();
	Det = 1;
	for (int i = 0; i < m.size(); ++i) {
		m[i][i] = 1;
	}
}

void Matrix::SetNull() {
	Det = 0;
	m.assign(m.size(), std::vector<double>(m.size(), 0));
}

double Matrix::det() {
	Matrix tmp(*this);
	for (int i = 0; i < (m.size() - 1); ++i) {
		if (abs(tmp.m[i][i]) < eps) {
			for (int j = i + 1; j < m.size(); ++j) {
				if (abs(tmp.m[j][i]) < eps) {
					tmp.m[i].swap(tmp.m[j]);
					break;
				}
				if (j == m.size() - 1)
					return Det = 0;
			}
		}

		for (int j = i + 1; j < m.size(); ++j) {
			for (int k = m.size() - 1; k >= i; --k) {
				tmp.m[j][k] -= tmp.m[j][i] * tmp.m[i][k] / tmp.m[i][i];
			}
		}
	}
	Det = 1.;
	for (int i = 0; i < m.size(); ++i) {
		Det *= tmp.m[i][i];
	}
	return Det;
}

std::vector<double> Matrix::operator*(std::vector<double>& b) {
	std::vector<double> res(m.size(), 0);
	for (int i = 0; i < b.size(); ++i) {
		for (int j = 0; j < b.size(); ++j) {
			res[i] += m[i][j] * b[j];
		}
	}
	return res;
}

Matrix Matrix::operator*(Matrix& b) {
	Matrix res(m.size());
	for (int i = 0; i < m.size(); ++i) {
		for (int j = 0; j < m.size(); ++j) {
			for (int k = 0; k < m.size(); ++k) {
				res.m[i][j] += m[i][k] * b.m[k][j];
			}
		}
	}
	return res;
}

Matrix& Matrix::operator=(Matrix& b) {
	Det = b.Det;
	m.assign(b.m.begin(), b.m.end());
	return *this;
}

std::pair<std::vector<double>, Matrix> Matrix::Gauss(std::vector<double> f) {
	Matrix t(*this);
	Matrix inverse(m.size());
	inverse.SetIdentity();
	//std::ofstream outf("output.txt");
	for (int i = 0; i < (m.size() - 1); ++i) {
		if (abs(m[i][i]) < eps) {
			for (int j = i + 1; j < m.size(); ++j) {
				if (abs(t.m[j][i]) < eps) {
					t.m[i].swap(t.m[j]);
					std::swap(f[i], f[j]);
					Matrix M(m.size());
					M.SetIdentity();
					M.m[i].swap(M.m[j]);
					Matrix tmp(M * inverse);
					inverse = tmp;
					break;
				}
			}
		}
		Matrix M(m.size());
		M.SetIdentity();
		for (int j = i + 1; j < m.size(); ++j)
			M.m[j][i] = -t.m[i][j] / t.m[i][i];
		Matrix tmp(M * inverse);
		inverse = tmp;
		/*for (int j = 0; j < 5; ++j){
			for (int k = 0; k < 5; ++k)
				outf << inverse.m[j][k] << ' ';
			outf << std::endl;
		}*/
		for (int j = i + 1; j < m.size(); ++j) {
			f[j] -= f[i] * t.m[j][i] / t.m[i][i];
			for (int k = m.size() - 1; k >= i; --k) {
				t.m[j][k] -= t.m[j][i] * t.m[i][k] / t.m[i][i];
			}
		}
		/*outf << std::endl << std::endl;
		for (int j = 0; j < 5; ++j){
			for (int k = 0; k < 5; ++k)
				outf << t.m[j][k] << ' ';
			outf << std::endl;
		}
		outf << std::endl << std::endl << std::endl << std::endl << std::endl;*/
	}
	Matrix M(m.size());
	M.SetNull();
	for (int i = m.size() - 1; i >= 0; --i) {
		M.m[i][i] = 1. / t.m[i][i];
	}
	for (int i = 0; i < m.size(); ++i) {
		for (int j = i + 1; j < m.size(); ++j) {
			double s = 0;
			for (int k = 0; k < j; ++k)
				s -= M.m[i][k] * t.m[k][j];
			M.m[i][j] = s * M.m[j][j];
		}
	}
	Matrix tmp(M * inverse);
	inverse = tmp;
	std::vector<double> res(m.size(), 0);
	res[m.size() - 1] = f[m.size() - 1] / t.m[m.size() - 1][m.size() - 1];
	for (int i = m.size() - 2; i >= 0; --i) {
		res[i] = f[i];
		for (int j = m.size() - 1; j > i; --j)
			res[i] -= t.m[i][j] * res[j];
		res[i] /= t.m[i][i];
	}
	return std::make_pair(res, inverse);
}

std::vector<double> Matrix::Jacobi(std::vector<double>& f) {
	std::vector <double> res(m.size(), 1.), rest(m.size());
	double change = 10.;
	while (change > eps) {
		for (int i = 0; i < m.size(); ++i) {
			rest[i] = f[i] / m[i][i];
			for (int j = 0; j < m.size(); ++j)
				rest[i] -= (i == j) ? 0. : res[j] * m[i][j] / m[i][i];
			change = std::min(change, std::abs(res[i] - rest[i]));
		}
		res.swap(rest);
	}
	return rest;
}

double Matrix::Norm() {
	double res = -1;
	for (int i = 0; i < m.size(); ++i){
		double s = 0.;
		for (int j = 0; j < m.size(); ++j){
			s += m[i][j] > 0 ? m[i][j] : -m[i][j];
		}
		res = std::max(s, res);
	}
	return res;
}


Matrix::~Matrix() {
}


Matrix Matrix::operator+(Matrix &b){
	Matrix res(m.size());
	for (int i = 0; i < m.size(); ++i) {
		for (int j = 0; j < m.size(); ++j) {
			res.m[i][j] += m[i][j] + b.m[i][j];
		}
	}
	return res;
}

Matrix Matrix::operator*(double b){
	Matrix res(*this);
	for (int i = 0; i < m.size(); ++i){
		for (int j = 0; j < m.size(); ++j)
			res.m[i][j] *= b;
	}
	return res;
}