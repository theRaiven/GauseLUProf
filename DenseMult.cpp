// DenseMult.cpp

#include "Header.h"

void MultDense( real** A, const int* xTrue, real* b, int n )
{
    for (int i = 0; i < n; i++)
    {
        realS sum{ 0.0 };
        for (int j = 0; j < n; j++)
        {
             sum += A[i][j] * xTrue[j];
        }
        b[i] = sum;
    }
}
void SubtractVectors(int n, const int* xTrue, const real* xFound)
{
    cout << endl << "x^* - x^k" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << std::fixed << std::setprecision(accuracy) << xTrue[i] - xFound[i] << endl;
    }
}
void WorlingWithDense()
{
    for (int k = 0; k <= 20; k++)
    {
        int n;
        real** A;
        int* xTrue;
    
        InputMatrix(A, xTrue, n);
        real original_di0 = A[0][0];

        std::cout << "\n=== Запуск #" << k << " ===" << std::endl;

        A[0][0] = original_di0 + pow(10.0, -k);
        //PrintMatrix(A, n);
        real* b = new real[n];
        MultDense(A, xTrue, b, n);


        real* x = new real[n];
        SolveGauss(A, b, x, n);

        SubtractVectors(n, xTrue, x);

        std::cout << std::endl << "x^k" << std::endl;
        Output(x, n);

        delete[] b;
        delete[] x;
        delete[] xTrue;

        for (int i = 0; i < n; i++)
        {
            delete[] A[i];
        }
    }
}