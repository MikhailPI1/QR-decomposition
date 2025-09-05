#ifndef QR_DECOMPOSITION_H
#define QR_DECOMPOSITION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <complex>

#define TOLERANCE 0.0001

using namespace std;

typedef vector<double> Vector;
typedef pair<complex<double>, complex<double>> ComplexPair;

class Matrix {
    vector<vector<double>> data;

public:
    Matrix(int n);
    Matrix(const Matrix& other);
    size_t size() const;
    void set(int i, int j, double x);
    vector<double> getRow(size_t i) const;
    double get(size_t i, size_t j) const;
    
    Matrix operator*(const Matrix& other) const;
    Matrix operator+(const Matrix& other) const;
    Matrix operator/(double scalar) const;
    Matrix operator*(double scalar) const;
    Matrix operator-(const Matrix& other) const;
    Matrix& operator=(const Matrix& other);
    
    void makeIdentity();
    void print() const;
    
    Matrix transpose() const;
    double norm() const;
};


Matrix outerProduct(const Vector& vec);
double innerProduct(const Vector& vec);
double sign(double x);

Matrix householder(Matrix& matrix, size_t k);
Matrix qrDecomposition(Matrix& matrix);
pair<vector<complex<double>>, Matrix> findEigenvalues(Matrix& matrix, int maxIter);

void matrixVectorMultiply(const Matrix& matrix, const Vector& vec, Vector& result);
void vectorScalarMultiply(const Vector& vec, complex<double> scalar, vector<complex<double>>& result);
void verifyEigenvectors(const vector<complex<double>>& eigenvalues, const Matrix& eigenvectors, const Matrix& original);

#endif