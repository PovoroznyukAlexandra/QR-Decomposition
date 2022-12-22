#include <iostream>
#include <ctime>
#include <cmath>
#include <chrono>
#include <vector>


class Matrix {
    std::vector<double> matrix;
    int n;
public:
    Matrix(int n, int m): n(n){
        matrix.resize(n * m);
    }

    double& operator()(int i, int j){
        return matrix[i * n + j];
    }

    double operator()(int i, int j) const{
        return  matrix[i * n + j];
    }
};

// Вращение матрицы
void Rotation(int n, int i, int j, double c, double s,  Matrix &A){
    for(int k = 0; k < n; k++){
        double Aik = c * A(i, k) - s * A(j, k);
        double Ajk = s * A(i, k) + c * A(j, k);
        A(i, k) = Aik;
        A(j, k) = Ajk;
    }
}

// Транспонирование матрицы Q
void Transpose(int n, Matrix &Q){
    double t;
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            t = Q(i, j);
            Q(i, j) = Q(j, i);
            Q(j, i) = t;
        }
    }
}

// Проверка того, что R - правая треугольная матрица
void Check_Right_Triang(int n, Matrix const &R){
    double eps = 1e-7;
    for (int i = 0; i < (n - 1); i++){
        for (int j = (i + 1); j < n; j++){
            if (R(j, i) > eps){
                std::cout << "Matrix R isn't right triangular!" << std::endl;
                exit(1);
            }
        }
    }
    std::cout << "Matrix R is right triangular!" << std::endl;
}

// Проверка унитарности Q
void Check_Unitary(int n, Matrix const &Q){
    Matrix C(n, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            C(i, j) = 0;
            for (int k = 0; k < n; k++) {
                C(i, j) += Q(i, k) * Q(j, k);
            }
        }
        C(i, i) -= 1;
    }
    double eps = 1e-0;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (C(i, j) > eps){
                std::cout << "Matrix Q isn't unitary!" << std::endl;
                return;
            }
        }
    }
    std::cout << "Matrix Q is unitary!" << std::endl;
}

// Подсчет ошибки вычисления
void Check_Accuracy(int n, Matrix const &A, Matrix const &Q, Matrix const &R){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double tmp = A(i, j);
            for (int k = 0; k < n; k++) {
                tmp -= Q(i, k) * R(k, j);
            }
            if (fabs(tmp) > 1e-10){
                std::cout << "QR decomposition isn't correct!" << std::endl;
                exit(1);
            }
        }
    }
    std::cout << "QR decomposition is correct!" << std::endl;
}


int main(){
    srand(time(nullptr));
    int n = 0;
    std::cout << "Enter N: ";
    std::cin >> n;
    std::cout << "╰( ͡° ͜ʖ ͡° )つ──☆*" << std::endl;

    // Исходная матрица
    Matrix A(n, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            A(i, j) = rand() % 100; // Каждый элемент случайному числу от 0 до 99
        }
    }

    // Унитарная матрица Q
    Matrix Q(n, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if(i == j){
                Q(i, j) = 1;
            }
            else{
                Q(i, j) = 0;
            }
        }
    }

    // Правая треугольная матрца R
    Matrix R(n, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            R(i, j) = A(i, j);
        }
    }

//-------------------------------------------------------------------------------------------------------------//
    // Начало замера времени
    auto begin = std::chrono::steady_clock::now();

    // Вращение матриц
    for (int i = 0; i < (n - 1); i++){
        for (int j = (i + 1); j < n; j++){
            double mod = sqrt(R(i, i) * R(i, i) + R(j, i) * R(j, i));
            double c = R(i, i) / mod, s = -R(j, i) / mod;
            Rotation(n, i, j, c, s, Q);
            Rotation(n, i, j, c, s, R);
        }
    }
    Transpose(n, Q);

    // Вывод времени работы программы
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - begin;
    std::cout << "Duration: " << elapsed_seconds.count() << " sec\n";
//-------------------------------------------------------------------------------------------------------------//

    // Проверка того, что R - правая треугольная матрица
    Check_Right_Triang(n, R);

    // Проверка унитарности Q
    Check_Unitary(n, Q);

    // Подсчет ошибки вычисления
    Check_Accuracy(n, A, Q, R);

    return 0;
}
