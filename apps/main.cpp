#define _USE_MATH_DEFINES
#include "qr_decomposition.h"
#include <iostream>
#include <string>

using namespace std;

int main() {
    try {
        int choice;
        cout << "Select the input method:" << endl;
        cout << "1 - Keyboard input" << endl;
        cout << "2 - Input from a file" << endl;
        cout << "Your choice: ";
        cin >> choice;
        
        int n, i, j;
        double temp;
        Matrix* matrix = nullptr;
        
        if (choice == 1) {
            cout << "Enter the size of the matrix n: ";
            cin >> n;
            matrix = new Matrix(n);
            
            cout << "Enter the elements of the matrix " << n << "x" << n << ":" << endl;
            for (i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    cout << "Element [" << i << "][" << j << "]: ";
                    cin >> temp;
                    matrix->set(i, j, temp);
                }
            }
            
            cout << endl << "Entered matrix:" << endl;
            matrix->print();
            
        } else if (choice == 2) {
            string filename;
            cout << "Input file name: ";
            cin >> filename;
            
            ifstream in(filename);
            if (!in.is_open()) {
                throw runtime_error("File is not open...");
            }
            
            in >> n;
            matrix = new Matrix(n);
            
            for (i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    in >> temp;
                    matrix->set(i, j, temp);
                }
            }
            
            in.close();
            
            cout << "Matrix from file:" << endl;
            matrix->print();
            
        } else {
            throw runtime_error("Wrong choice");
        }
        
        // Create copy for calculations
        Matrix copy_matrix = *matrix;
        
        // Find eigenvalues and eigenvectors
        pair<vector<complex<double>>, Matrix> result = findEigenvalues(copy_matrix, n * 50);
        
        // Display results
        cout << "\n=== RESULTS ===" << endl;
        for (int i = 0; i < n; i++) {
            cout << "Eigenvalue " << i+1 << ": " << result.first[i] << endl;
            cout << "Eigenvector " << i+1 << ": ";
            Vector eigenvector = result.second.getRow(i);
            for (int j = 0; j < n; j++) {
                cout << eigenvector[j] << " ";
            }
            cout << endl << endl;
        }
        
        // Verification
        cout << "=== VERIFICATION ===" << endl;
        verifyEigenvectors(result.first, result.second, *matrix);
        
        // Clean up
        delete matrix;
        
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    } catch (const exception& e) {
        cerr << "Unexpected error: " << e.what() << endl;
        return 1;
    }

    return 0; 
}