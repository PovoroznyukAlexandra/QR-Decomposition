#include <iostream>
#include <chrono>
#include <vector>

using namespace std;

extern "C" void dgeqrf_(int* m, int* n, double* A, int* LDA, double* tau, double* work, int* lwork, int* info);

int main(){
    int n = 0;
    cout << "Enter N: ";
    cin >> n;
    cout << "╰( ͡° ͜ʖ ͡° )つ──☆*" << endl;

    int LDA = n;
    int lwork = -1;
    int info = 0;
    auto* tau = new double[n];
    auto* l_work = new double[n];

    // Исходная матрица
    auto* A = new double[n*n];
    for(int i = 0; i < n * n; i++) {
        A[i] = int(rand()) % 100;
    }

    dgeqrf_(&n, &n, A, &LDA, tau, l_work, &lwork, &info);
    lwork = l_work[0];

    auto* work = new double[lwork];
    auto start = chrono::steady_clock::now();
    dgeqrf_(&n, &n, A, &LDA, tau, work, &lwork, &info);

    auto end = chrono::steady_clock::now();
    cout << "Duration: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << "\n";
    delete[] tau;
    delete[] work;
    delete[] l_work;
//-------------------------------------------------------------------------------------------------------------//
    return 0;
}
