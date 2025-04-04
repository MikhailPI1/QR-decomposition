#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <complex>

#define tolerance 0.001 //Наша приближение

using namespace std;

ifstream in("matrix.txt");

//Общий класс матриц размером n*n
class Matrix {
    vector<vector<double>> data;

    public:

    Matrix(int n) :
    data(n, vector<double>(n))
    {}
    
    Matrix(const Matrix& other, int n) : data(n) {
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
            data[i][j] = other.data[i][j];
            }
        }
    }


    size_t size() const{
        //Вернет размер матрицы
        return data.size();
    }

    void pull(int i, int j, double x) {
        //Функция заполнения матрицы 
        data[i][j] = x;
    }

    Matrix operator*(Matrix& other) {

        //Перегрузка оператора для перемножения матриц.
        // Соотвественно вернет матрицы n*n; При попытки перемножить матрицы разных размерностей выдаст ошибку
        try {

            size_t n = data.size(), i, j, k;
            if (other.size() != n) {
                throw runtime_error("Матрицы имеют разный размер");
            }

            Matrix new_matrix(n);

            for (i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    for (k = 0; k < n; ++k) {
                        new_matrix.data[i][j] += data[i][k] * other.data[k][j];
                    }
                }
            }

            return new_matrix;

        } catch (const runtime_error& e) {
            cerr << e.what() << endl;
            return 0;
        } 
    }
    vector <double> get(size_t i) {
        return data[i];
    }

    double get(size_t i, size_t j){
        return data[i][j];
    }

    Matrix operator+(Matrix& other) {

        //Перегрузка оператора для сложения матриц.
        // Соотвественно вернет матрицы n*n; При попытки сложения матрицы разных размерностей выдаст ошибку
        try {
            size_t n = data.size(), i, j;
            if (other.size() != n) {
                throw runtime_error("Матрицы имеют разный размер");
            }

            Matrix new_matrix(n);

            for (i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    new_matrix.data[i][j] = data[i][j] - other.data[i][j];
                }
            }

            return new_matrix;

        } catch (const runtime_error& e) {
            cerr << e.what() << endl;
        }
    }

    Matrix operator/(double scalar) {
        //Перегрузка оператора для деления матрицы на скаляр.
        // При попытке деления на 0 выста ошибку
        try {
            if (scalar == 0) {
                throw runtime_error("Попытка делить на ноль");
            }
            size_t i, j, n = data.size();

            for (i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    data[i][j] = data[i][j] / scalar;
                }
            }
            return *this;
        } catch (const runtime_error& e){
            cerr << e.what() << endl;
        }
    }

    Matrix operator*(double scalar) {
        //Перегрузка оператора для умножения матрицы на скаляр.

        size_t i, j, n = data.size();

        for (i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                data[i][j] = data[i][j] * scalar;
            }
        }
        return *this;
    }

    Matrix operator-(Matrix& other) {
        //Перегрузка оператора для нахождения разности матриц;
        // Соотвественно вернет матрицы n*n; При попытки вычесть разницу матриц разных размерностей выдаст ошибку

        try {
            size_t n = data.size(), i, j;
            if (other.size() != n) {
                throw runtime_error("Матрицы имеют разный размер");
            }

            Matrix new_matrix(n);

            for (i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    new_matrix.data[i][j] = data[i][j] - other.data[i][j];
                }
            }

            return new_matrix;

        } catch (const runtime_error& e) {
            cerr << e.what() << endl;
            return 0;
        }
    }

    Matrix& operator=(const Matrix& other)  {
        //Оператор присваивания.

        if (this == &other) {
            return *this;
          }
      
          Matrix temp(other);
      
          swap(data, temp.data);

          return *this;
    }

    void Singular() {
        //Превращает матрицу в единичную.

        size_t n = data.size(), i, j;
        for (i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                if (i == j) {
                    data[i][j] = 1;
                } else {
                    data[i][j] = 0;
                }
            }
        } 
    }
    
    void view() {

        cout << endl;
        //Выведет текущий вид матрицы.
        size_t i, j, n = data.size();
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }
};

typedef vector<double> Vector; //Общее обозначение вектора
typedef pair<complex<double>, complex<double>> complex_pair; //Пара комплексных чисел.


Matrix vectorTranspose(const Vector& vec) {
    //Функция перемножения векторов в порядке v * vT. Выдаст матрицу.
    size_t n = vec.size(), i, j;

    Matrix resultMatrix(n);

    for ( i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            resultMatrix.pull(i, j, vec[i] * vec[j]);
        }
    }

    return resultMatrix;
}

double transposeVector(const Vector& vec) {

    //Функция перемножения векторов в порядке vT * v. Выдаст скаляр.
    int n = vec.size();
    double scalar = 0;

    for (int i = 0; i < n; ++i) {
        scalar += vec[i] * vec[i];
    }

    return scalar;
}

double sign(double x) {

    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;

    return 0.0;
}
double sign(int x) {

    if (x > 0) return 1;
    if (x < 0) return -1;
    
    return 0.0;
}
complex_pair diskreminant(Matrix matrix, size_t i){

    complex_pair self_value;

    double a = matrix.get(i, i),
           b = matrix.get(i + 1, i + 1),
           c = matrix.get(i + 1, i),
           d = matrix.get(i, i + 1);

    double disk = (a + b) * (a + b) - 4 * (a * b - c * d);

    if (disk >= 0) {
        // Вещественные корни
        double root1 = ( (a + b) + sqrt(disk)) / 2.0;
        double root2 = ((a + b) - sqrt(disk)) / 2.0;

        return {root1, root2};
    } else {
        // Комплексные корни
        complex<double> realPart = (a + b) / 2.0;
        complex<double> imaginaryPart = sqrt(complex<double>(disk)) / (2.0);
        complex<double> root1 = realPart + imaginaryPart;
        complex<double> root2 = realPart - imaginaryPart;
        return {root1, root2};
    }
}


Matrix haus(Matrix& matrix, size_t k) {

    //Данная функция находит матрицу Хаусхолдера.

    size_t n = matrix.size(), i;
    Matrix E(n);
    Vector vec(n);
    double elem = 0;

    E.Singular();

    for (i = 0; i < k; i++) {
        vec[i] = 0;
    }

    for(i = k + 1; i < n; i++) {
        vec[i] = matrix.get(i, k);
        elem += matrix.get(i,k) * matrix.get(i, k);
    }
    cout << endl;

    elem += matrix.get(k, k) * matrix.get(k, k);
    elem = sqrt(elem);
    vec[k] = matrix.get(k, k) + sign(matrix.get(k, k)) * elem;

    Matrix transVec(n);
    transVec = (vectorTranspose(vec) * 2) / transposeVector(vec);

    
    E = E  - transVec;
    
    matrix = E * matrix;
    return E;
}

Matrix QR (Matrix& matrix) {
    
    //Данная функция выполняет QR-разложения матрицы и возваращает ортогональную матрицу Q;
    //Так же происхожит обратное перемножение RQ что бы найти матрицу A_i+1

    size_t i, n = matrix.size();
    Matrix Q = haus(matrix, 0), temp(n);
    
    for(i = 1; i  < n - 1; i++) {
        temp = haus(matrix, i);
        Q = Q * temp;
    }
    
    matrix = matrix * Q;
    matrix.view();
    return Q;
}


pair<vector<complex<double>>, Matrix> selfValue(Matrix& matrix, int maxIter) {

    //Данная функция приводит исходную матрицу  блочнотреугольному ввиду(при наличии комплексных корней) или к верхнетреугольному виду.
    //А так же возвращает накопленную матрицу Q, столбцы которой являются собственными векторами.

    size_t i, j, n = matrix.size();
    int counter = 0;

    Matrix Q(n);
    Q.Singular();
    //Сохраняя ортогональную матрицу из прошлого разложения

    Matrix temp(n);

    vector<complex<double>> selfvalue; //Вектор возвращаемых значений

    for(int iter = 0; iter < maxIter; iter++) { //Данный цикл служит для приведения матриц к верхне-треугольному виду.

        temp = QR(matrix);
        Q = Q * temp; //Накапливаем Q;

        bool converged = true;
        for (i = 0; i < n - 1; ++i) {
            for (j = i + 1; j < n; ++j) {
                if (abs(matrix.get(j, i)) > tolerance) {
                    converged = false;
                    break;
                }
            }
            if (!converged) break;
        }

        if (converged) {
            break; 
        }
    }

    for(i = 0; i < n; i++) { //Данный цикл служит для вычленения собственных значения и нахождения комплексных пар.
        if (n - i == 1) {
            complex<double> value;
            value = matrix.get(i,i);
            selfvalue.push_back(value);
            break;
        }
        if(matrix.get(i + 1, i) >= 0.0001) {
            counter++;
            if (counter == 1) {
                counter = 0;
                complex_pair pair = diskreminant(matrix, i);
                selfvalue.push_back(pair.first);
                selfvalue.push_back(pair.second);
            }
        } else {
            complex<double> value;
            value = matrix.get(i,i);
            selfvalue.push_back(value);
            counter = 0;
        }
    }
    pair<vector<complex<double>>, Matrix> return_value = make_pair(selfvalue, Q);
    return return_value;
}

void matrix_on_vector(Matrix matrix, vector<double> vector) {
    size_t n = matrix.size();
    //Функция перемножения матрицы на вектор
    Vector result(n);
    for (size_t i = 0; i < n; ++i) {
        result[i] = 0.0;
      for (size_t j = 0; j < vector.size(); ++j) {
        result[i] += matrix.get(i, j) * vector[j];
      }
    }

    for (size_t i = 0; i < matrix.size(); i++) {
        cout << result[i] << " ";
    }
    cout << endl;
}
void vector_on_scalar(vector<double> vector, complex<double> scalar) {

    //Функция перемножения вектора на скаляр;
    complex<double> temp;

    for(size_t i = 0; i < vector.size(); i++) {
        temp = vector[i];
        cout << temp * scalar << " ";
    }
    cout << endl;
}

void cheak(vector<complex<double>> selfvalue, Matrix Q, Matrix matrix) {
    cout << endl;
    for (size_t i = 0; i < matrix.size(); i++) {
        cout << "Произведение вектора " << i << " собственное значение " << endl;
        vector_on_scalar(Q.get(i), selfvalue[i]);

        cout << "Произведение матрицы на собственный вектор " << i <<endl;
        matrix_on_vector(matrix, Q.get(i));
        cout << endl;
    }
}

int main() {

    try {
        if (!in.is_open()) {
            throw runtime_error("Файл не был открыт...");
        }

        int n, i, j, temp;

        in >> n;
        Matrix matrix(n);
        

        for (i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                in >> temp;

                matrix.pull(i, j, temp);
            }
        }
        
        Matrix copy_matrix = matrix;

        pair<vector<complex<double>>, Matrix> pair = selfValue(copy_matrix, n * n);
        cout << endl;
        for (int i = 0; i < n; i++) {
            cout << "Собственное значение " <<pair.first[i]  << endl;
            cout << "Собственный вектор ";
            for (int j = 0; j < n; j++) {
                cout << pair.second.get(i, j) << " ";
            }
            cout << endl;
        }

        cheak(pair.first, pair.second, matrix);

    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
    }

   return 0; 
}