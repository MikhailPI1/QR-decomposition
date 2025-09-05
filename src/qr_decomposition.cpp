#include "qr_decomposition.h"

extern ifstream in;

Matrix::Matrix(int n) : data(n, vector<double>(n, 0.0)) {}

Matrix::Matrix(const Matrix& other) : data(other.data) {}

size_t Matrix::size() const {
    return data.size();
}

void Matrix::set(int i, int j, double x) {
    data[i][j] = x;
}

vector<double> Matrix::getRow(size_t i) const {
    return data[i];
}

double Matrix::get(size_t i, size_t j) const {
    return data[i][j];
}

Matrix Matrix::operator*(const Matrix& other) const {
    size_t n = data.size();
    if (other.size() != n) {
        throw runtime_error("Matrices have different sizes");
    }

    Matrix result(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < n; k++) {
                result.data[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }
    return result;
}

Matrix Matrix::operator+(const Matrix& other) const {
    size_t n = data.size();
    if (other.size() != n) {
        throw runtime_error("Matrices have different sizes");
    }

    Matrix result(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            result.data[i][j] = data[i][j] + other.data[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator/(double scalar) const {
    if (abs(scalar) < 1e-10) {
        throw runtime_error("Division by zero");
    }

    Matrix result(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data.size(); j++) {
            result.data[i][j] = data[i][j] / scalar;
        }
    }
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data.size(); j++) {
            result.data[i][j] = data[i][j] * scalar;
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    size_t n = data.size();
    if (other.size() != n) {
        throw runtime_error("Matrices have different sizes");
    }

    Matrix result(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            result.data[i][j] = data[i][j] - other.data[i][j];
        }
    }
    return result;
}

Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        data = other.data;
    }
    return *this;
}

void Matrix::makeIdentity() {
    size_t n = data.size();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            data[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

void Matrix::print() const {
    cout << endl;
    size_t n = data.size();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

Matrix Matrix::transpose() const {
    size_t n = data.size();
    Matrix result(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            result.data[j][i] = data[i][j];
        }
    }
    return result;
}

double Matrix::norm() const {
    double sum = 0.0;
    for (const auto& row : data) {
        for (double val : row) {
            sum += val * val;
        }
    }
    return sqrt(sum);
}

Matrix outerProduct(const Vector& vec) {
    size_t n = vec.size();
    Matrix result(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            result.set(i, j, vec[i] * vec[j]);
        }
    }
    return result;
}

double innerProduct(const Vector& vec) {
    double result = 0.0;
    for (double val : vec) {
        result += val * val;
    }
    return result;
}

double sign(double x) {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

ComplexPair solveQuadratic(double a, double b, double c, double d) {
    double trace = a + d;
    double determinant = a * d - b * c;
    double discriminant = trace * trace - 4 * determinant;
    
    if (discriminant >= 0) {
        double root1 = (trace + sqrt(discriminant)) / 2.0;
        double root2 = (trace - sqrt(discriminant)) / 2.0;
        return {root1, root2};
    } else {
        complex<double> realPart = trace / 2.0;
        complex<double> imaginaryPart = sqrt(-discriminant) / 2.0;
        return {realPart + imaginaryPart, realPart - imaginaryPart};
    }
}

Matrix householder(Matrix& matrix, size_t k) {
    size_t n = matrix.size();
    Vector x(n, 0.0);
    
    for (size_t i = k; i < n; i++) {
        x[i] = matrix.get(i, k);
    }
    
    double norm_x = sqrt(innerProduct(x));
    x[k] += sign(x[k]) * norm_x;
    
    double beta = innerProduct(x);
    if (abs(beta) < 1e-10) {
        Matrix I(n);
        I.makeIdentity();
        return I;
    }
    
    Matrix H = outerProduct(x) * (2.0 / beta);
    Matrix I(n);
    I.makeIdentity();
    H = I - H;
    
    matrix = H * matrix;
    return H;
}

Matrix qrDecomposition(Matrix& matrix) {
    size_t n = matrix.size();
    Matrix Q(n);
    Q.makeIdentity();
    
    for (size_t k = 0; k < n - 1; k++) {
        Matrix H = householder(matrix, k);
        Q = Q * H;
    }
    
    return Q;
}

pair<vector<complex<double>>, Matrix> findEigenvalues(Matrix& matrix, int maxIter) {
    size_t n = matrix.size();
    Matrix Q_total(n);
    Q_total.makeIdentity();
    
    vector<complex<double>> eigenvalues;
    
    for (int iter = 0; iter < maxIter; iter++) {
        Matrix Q = qrDecomposition(matrix);
        Q_total = Q_total * Q;
        matrix = matrix * Q;
        
        bool converged = true;
        for (size_t i = 1; i < n; i++) {
            for (size_t j = 0; j < i; j++) {
                if (abs(matrix.get(i, j)) > TOLERANCE) {
                    converged = false;
                    break;
                }
            }
            if (!converged) break;
        }
        
        if (converged) break;
    }
    
    size_t i = 0;
    while (i < n) {
        if (i == n - 1 || abs(matrix.get(i + 1, i)) < TOLERANCE) {
            eigenvalues.push_back(matrix.get(i, i));
            i++;
        } else {
            ComplexPair roots = solveQuadratic(
                matrix.get(i, i), matrix.get(i, i + 1),
                matrix.get(i + 1, i), matrix.get(i + 1, i + 1)
            );
            eigenvalues.push_back(roots.first);
            eigenvalues.push_back(roots.second);
            i += 2;
        }
    }
    
    return make_pair(eigenvalues, Q_total);
}

void matrixVectorMultiply(const Matrix& matrix, const Vector& vec, Vector& result) {
    size_t n = matrix.size();
    result.resize(n, 0.0);
    
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            result[i] += matrix.get(i, j) * vec[j];
        }
    }
}

void vectorScalarMultiply(const Vector& vec, complex<double> scalar, vector<complex<double>>& result) {
    result.resize(vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
        result[i] = vec[i] * scalar;
    }
}

void verifyEigenvectors(const vector<complex<double>>& eigenvalues, const Matrix& eigenvectors, const Matrix& original) {
    cout << "\n=== Eigenvector Verification ===\n";
    size_t n = original.size();
    
    for (size_t i = 0; i < n; i++) {
        cout << "\nEigenvalue " << i + 1 << ": " << eigenvalues[i] << endl;
        
        Vector eigenvector = eigenvectors.getRow(i);
        
        Vector Av;
        matrixVectorMultiply(original, eigenvector, Av);
        
        vector<complex<double>> lambdaV;
        vectorScalarMultiply(eigenvector, eigenvalues[i], lambdaV);
        
        cout << "A * v: ";
        for (double val : Av) cout << val << " ";
        
        cout << "\nλ * v: ";
        for (const auto& val : lambdaV) cout << val << " ";
        
        double maxError = 0.0;
        for (size_t j = 0; j < n; j++) {
            double error = abs(Av[j] - real(lambdaV[j]));
            maxError = max(maxError, error);
        }
        
        cout << "\nMax error: " << maxError << endl;
        cout << "Status: " << (maxError < TOLERANCE ? "PASS" : "FAIL") << endl;
    }
}