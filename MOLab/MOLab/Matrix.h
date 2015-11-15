#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
class Matrix {
public:
	Matrix();
	Matrix(int n);
	void SetIdentity();
	void SetNull();
	Matrix(const Matrix& orig);
	virtual ~Matrix();
	Matrix& operator=(Matrix &b);
	std::vector< std::vector<double> > m; \
	std::vector<double> operator*(std::vector<double> &b);
	std::pair< std::vector<double>, Matrix> Gauss(std::vector<double> f);
	std::vector <double> Jacobi(std::vector<double> &f);
	Matrix operator*(Matrix &b);
	Matrix operator*(double b);
	Matrix operator+(Matrix &b);
	double Norm();
	double det();
	double Det;
private:

};
