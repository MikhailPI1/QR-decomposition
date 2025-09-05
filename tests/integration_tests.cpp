#include "qr_decomposition.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

// Integration test for complete QR decomposition workflow
void testQRDecompositionWorkflow() {
    cout << "=== Integration Test: QR Decomposition Workflow ===" << endl;
    

    Matrix A(3);
    A.pull(0, 0, 12.0);
    A.pull(0, 1, -51.0);
    A.pull(0, 2, 4.0);
    A.pull(1, 0, 6.0);
    A.pull(1, 1, 167.0);
    A.pull(1, 2, -68.0);
    A.pull(2, 0, -4.0);
    A.pull(2, 1, 24.0);
    A.pull(2, 2, -41.0);
    
    cout << "Original matrix A:" << endl;
    A.view();
    

    Matrix original = A;
    Matrix Q = QR(A);
    Matrix R = A; 
    
    cout << "Q matrix:" << endl;
    Q.view();
    
    cout << "R matrix:" << endl;
    R.view();
    

    Matrix reconstructed = Q * R;
    cout << "Reconstructed matrix Q * R:" << endl;
    reconstructed.view();
    

    double max_error = 0.0;
    for (size_t i = 0; i < original.size(); ++i) {
        for (size_t j = 0; j < original.size(); ++j) {
            double error = abs(original.get(i, j) - reconstructed.get(i, j));
            max_error = max(max_error, error);
        }
    }
    
    cout << "Maximum reconstruction error: " << max_error << endl;
    cout << "Test " << (max_error < 0.001 ? "PASSED" : "FAILED") << endl << endl;
}

// Integration test for eigenvalue computation
void testEigenvalueComputation() {
    cout << "=== Integration Test: Eigenvalue Computation ===" << endl;
    

    Matrix A(3);
    A.pull(0, 0, 2.0);
    A.pull(0, 1, -1.0);
    A.pull(0, 2, 0.0);
    A.pull(1, 0, -1.0);
    A.pull(1, 1, 2.0);
    A.pull(1, 2, -1.0);
    A.pull(2, 0, 0.0);
    A.pull(2, 1, -1.0);
    A.pull(2, 2, 2.0);
    
    cout << "Symmetric matrix A:" << endl;
    A.view();
    

    auto result = selfValue(A, 50);
    vector<complex<double>> eigenvalues = result.first;
    Matrix eigenvectors = result.second;
    
    cout << "Computed eigenvalues:" << endl;
    for (size_t i = 0; i < eigenvalues.size(); ++i) {
        cout << "λ" << i << " = " << eigenvalues[i] << endl;
    }
    

    vector<double> expected = {2.0, 2.0 + sqrt(2.0), 2.0 - sqrt(2.0)};
    sort(expected.begin(), expected.end());
    
    vector<double> computed_real;
    for (const auto& eval : eigenvalues) {
        if (abs(eval.imag()) < tolerance) {
            computed_real.push_back(eval.real());
        }
    }
    sort(computed_real.begin(), computed_real.end());
    
    cout << "Expected eigenvalues: ";
    for (double val : expected) cout << val << " ";
    cout << endl;
    
    cout << "Computed real eigenvalues: ";
    for (double val : computed_real) cout << val << " ";
    cout << endl;
    

    bool test_passed = true;
    for (size_t i = 0; i < expected.size(); ++i) {
        if (abs(expected[i] - computed_real[i]) > 0.1) {
            test_passed = false;
            break;
        }
    }
    
    cout << "Eigenvalue test " << (test_passed ? "PASSED" : "FAILED") << endl << endl;
}

// Test with known eigenvalues (diagonal matrix)
void testDiagonalMatrix() {
    cout << "=== Integration Test: Diagonal Matrix ===" << endl;
    
    Matrix A(3);
    A.pull(0, 0, 1.0);
    A.pull(0, 1, 0.0);
    A.pull(0, 2, 0.0);
    A.pull(1, 0, 0.0);
    A.pull(1, 1, 2.0);
    A.pull(1, 2, 0.0);
    A.pull(2, 0, 0.0);
    A.pull(2, 1, 0.0);
    A.pull(2, 2, 3.0);
    
    cout << "Diagonal matrix A:" << endl;
    A.view();
    
    auto result = selfValue(A, 10);
    vector<complex<double>> eigenvalues = result.first;
    
    cout << "Computed eigenvalues:" << endl;
    for (const auto& eval : eigenvalues) {
        cout << eval << endl;
    }
    

    vector<double> expected = {1.0, 2.0, 3.0};
    vector<double> computed;
    
    for (const auto& eval : eigenvalues) {
        if (abs(eval.imag()) < tolerance) {
            computed.push_back(eval.real());
        }
    }
    
    sort(computed.begin(), computed.end());
    
    bool test_passed = (computed.size() == 3);
    for (size_t i = 0; i < computed.size() && test_passed; ++i) {
        test_passed = (abs(computed[i] - expected[i]) < tolerance);
    }
    
    cout << "Diagonal matrix test " << (test_passed ? "PASSED" : "FAILED") << endl << endl;
}

int main() {
    cout << "Running Integration Tests for QR Decomposition" << endl;
    cout << "=============================================" << endl << endl;
    
    testQRDecompositionWorkflow();
    testEigenvalueComputation();
    testDiagonalMatrix();
    
    cout << "All integration tests completed." << endl;
    return 0;
}