#include "qr_decomposition.h"
#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

// Test Matrix construction and basic operations
void testConstructorAndSize() {
    cout << "Testing constructor and size..." << endl;
    Matrix m(5);
    assert(m.size() == 5);
    cout << "✓ PASSED" << endl;
}

void testGetAndPull() {
    cout << "Testing get and pull..." << endl;
    Matrix m(2);
    m.pull(0, 0, 5.0);
    m.pull(1, 1, 10.0);
    
    assert(abs(m.get(0, 0) - 5.0) < tolerance);
    assert(abs(m.get(1, 1) - 10.0) < tolerance);
    cout << "✓ PASSED" << endl;
}

void testMatrixMultiplication() {
    cout << "Testing matrix multiplication..." << endl;
    Matrix A(2);
    A.pull(0, 0, 1.0);
    A.pull(0, 1, 2.0);
    A.pull(1, 0, 3.0);
    A.pull(1, 1, 4.0);
    
    Matrix B(2);
    B.pull(0, 0, 5.0);
    B.pull(0, 1, 6.0);
    B.pull(1, 0, 7.0);
    B.pull(1, 1, 8.0);
    
    Matrix C = A * B;
    
    assert(abs(C.get(0, 0) - 19.0) < tolerance);
    assert(abs(C.get(0, 1) - 22.0) < tolerance);
    assert(abs(C.get(1, 0) - 43.0) < tolerance);
    assert(abs(C.get(1, 1) - 50.0) < tolerance);
    cout << "✓ PASSED" << endl;
}

void testMatrixAddition() {
    cout << "Testing matrix addition..." << endl;
    Matrix A(2);
    A.pull(0, 0, 1.0);
    A.pull(0, 1, 2.0);
    A.pull(1, 0, 3.0);
    A.pull(1, 1, 4.0);
    
    Matrix B(2);
    B.pull(0, 0, 5.0);
    B.pull(0, 1, 6.0);
    B.pull(1, 0, 7.0);
    B.pull(1, 1, 8.0);
    
    Matrix C = A + B;
    
    assert(abs(C.get(0, 0) - 6.0) < tolerance);
    assert(abs(C.get(0, 1) - 8.0) < tolerance);
    assert(abs(C.get(1, 0) - 10.0) < tolerance);
    assert(abs(C.get(1, 1) - 12.0) < tolerance);
    cout << "✓ PASSED" << endl;
}

void testMatrixSubtraction() {
    cout << "Testing matrix subtraction..." << endl;
    Matrix A(2);
    A.pull(0, 0, 5.0);
    A.pull(0, 1, 6.0);
    A.pull(1, 0, 7.0);
    A.pull(1, 1, 8.0);
    
    Matrix B(2);
    B.pull(0, 0, 1.0);
    B.pull(0, 1, 2.0);
    B.pull(1, 0, 3.0);
    B.pull(1, 1, 4.0);
    
    Matrix C = A - B;
    
    assert(abs(C.get(0, 0) - 4.0) < tolerance);
    assert(abs(C.get(0, 1) - 4.0) < tolerance);
    assert(abs(C.get(1, 0) - 4.0) < tolerance);
    assert(abs(C.get(1, 1) - 4.0) < tolerance);
    cout << "✓ PASSED" << endl;
}

void testScalarMultiplication() {
    cout << "Testing scalar multiplication..." << endl;
    Matrix A(2);
    A.pull(0, 0, 1.0);
    A.pull(0, 1, 2.0);
    A.pull(1, 0, 3.0);
    A.pull(1, 1, 4.0);
    
    Matrix B = A * 2.0;
    
    assert(abs(B.get(0, 0) - 2.0) < tolerance);
    assert(abs(B.get(0, 1) - 4.0) < tolerance);
    assert(abs(B.get(1, 0) - 6.0) < tolerance);
    assert(abs(B.get(1, 1) - 8.0) < tolerance);
    cout << "✓ PASSED" << endl;
}

void testScalarDivision() {
    cout << "Testing scalar division..." << endl;
    Matrix A(2);
    A.pull(0, 0, 2.0);
    A.pull(0, 1, 4.0);
    A.pull(1, 0, 6.0);
    A.pull(1, 1, 8.0);
    
    Matrix B = A / 2.0;
    
    assert(abs(B.get(0, 0) - 1.0) < tolerance);
    assert(abs(B.get(0, 1) - 2.0) < tolerance);
    assert(abs(B.get(1, 0) - 3.0) < tolerance);
    assert(abs(B.get(1, 1) - 4.0) < tolerance);
    cout << "✓ PASSED" << endl;
}

void testSingularMatrix() {
    cout << "Testing singular matrix..." << endl;
    Matrix A(3);
    A.Singular();
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) {
                assert(abs(A.get(i, j) - 1.0) < tolerance);
            } else {
                assert(abs(A.get(i, j) - 0.0) < tolerance);
            }
        }
    }
    cout << "✓ PASSED" << endl;
}

// Test utility functions
void testSignFunction() {
    cout << "Testing sign function..." << endl;
    assert(abs(sign(5.0) - 1.0) < tolerance);
    assert(abs(sign(-3.0) - (-1.0)) < tolerance);
    assert(abs(sign(0.0) - 0.0) < tolerance);
    assert(abs(sign(10) - 1.0) < tolerance);
    assert(abs(sign(-7) - (-1.0)) < tolerance);
    cout << "✓ PASSED" << endl;
}

void testVectorTranspose() {
    cout << "Testing vector transpose..." << endl;
    Vector vec = {1.0, 2.0, 3.0};
    Matrix result = vectorTranspose(vec);
    
    assert(abs(result.get(0, 0) - 1.0) < tolerance);
    assert(abs(result.get(0, 1) - 2.0) < tolerance);
    assert(abs(result.get(0, 2) - 3.0) < tolerance);
    assert(abs(result.get(1, 0) - 2.0) < tolerance);
    assert(abs(result.get(1, 1) - 4.0) < tolerance);
    assert(abs(result.get(1, 2) - 6.0) < tolerance);
    assert(abs(result.get(2, 0) - 3.0) < tolerance);
    assert(abs(result.get(2, 1) - 6.0) < tolerance);
    assert(abs(result.get(2, 2) - 9.0) < tolerance);
    cout << "✓ PASSED" << endl;
}

void testTransposeVector() {
    cout << "Testing transpose vector..." << endl;
    Vector vec = {1.0, 2.0, 3.0};
    double result = transposeVector(vec);
    
    assert(abs(result - 14.0) < tolerance);
    cout << "✓ PASSED" << endl;
}

// Test Householder transformation
void testHouseholderTransformation() {
    cout << "Testing Householder transformation..." << endl;
    Matrix A(3);
    A.pull(0, 0, 1.0);
    A.pull(0, 1, 2.0);
    A.pull(0, 2, 3.0);
    A.pull(1, 0, 4.0);
    A.pull(1, 1, 5.0);
    A.pull(1, 2, 6.0);
    A.pull(2, 0, 7.0);
    A.pull(2, 1, 8.0);
    A.pull(2, 2, 9.0);
    
    Matrix original = A;
    Matrix H = haus(A, 0);
    
    // Check that transformation occurred
    assert(abs(A.get(1, 0) - original.get(1, 0)) > tolerance);
    assert(abs(A.get(2, 0) - original.get(2, 0)) > tolerance);
    cout << "✓ PASSED" << endl;
}

// Test QR decomposition
void testQRDecomposition() {
    cout << "Testing QR decomposition..." << endl;
    Matrix A(2);
    A.pull(0, 0, 3.0);
    A.pull(0, 1, 1.0);
    A.pull(1, 0, 1.0);
    A.pull(1, 1, 2.0);
    
    Matrix original = A;
    Matrix Q = QR(A);
    
    // Check that R is upper triangular
    assert(abs(A.get(1, 0)) < tolerance);
    cout << "✓ PASSED" << endl;
}

// Test eigenvalue computation for simple 2x2 matrix
void testEigenvalueComputation() {
    cout << "Testing eigenvalue computation..." << endl;
    Matrix A(2);
    A.pull(0, 0, 4.0);
    A.pull(0, 1, 1.0);
    A.pull(1, 0, 1.0);
    A.pull(1, 1, 3.0);
    
    auto result = selfValue(A, 100);
    vector<complex<double>> eigenvalues = result.first;
    
    // Expected eigenvalues: (7 ± √5)/2 ≈ 4.618 and 2.382
    assert(eigenvalues.size() == 2);
    
    bool found1 = false, found2 = false;
    for (const auto& eval : eigenvalues) {
        if (abs(eval - complex<double>(4.618, 0.0)) < 0.1) found1 = true;
        if (abs(eval - complex<double>(2.382, 0.0)) < 0.1) found2 = true;
    }
    
    assert(found1);
    assert(found2);
    cout << "✓ PASSED" << endl;

}

int main() {
    cout << "Running unit tests..." << endl;
    cout << "=====================" << endl;
    
    try {
        testConstructorAndSize();
        testGetAndPull();
        testMatrixMultiplication();
        testMatrixAddition();
        testMatrixSubtraction();
        testScalarMultiplication();
        testScalarDivision();
        testSingularMatrix();
        testSignFunction();
        testVectorTranspose();
        testTransposeVector();
        testHouseholderTransformation();
        testQRDecomposition();
        testEigenvalueComputation();
        
        cout << endl << "=====================" << endl;
        cout << "All tests PASSED! ✓" << endl;
        return 0;
        
    } catch (const exception& e) {
        cout << endl << "=====================" << endl;
        cout << "Test FAILED: " << e.what() << endl;
        return 1;
    }
}