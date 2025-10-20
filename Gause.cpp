// Gause.cpp

#include "Header.h"
// метод решения плотной матрицы по Гауссу
void SolveGauss(real** A, real* b, real*& x, int& n)
{
    try
    {
        const real EPS = relativeEPS<real>();

        for (int i = 0; i < n; i++)
        {
            int i0 = i;
            real maxElem = fabs(A[i][i]);

            for (int j = i + 1; j < n; j++)
            {
                if (fabs(A[j][i]) > maxElem) // 1 модуль, 1 сравнение 
                {
                    maxElem = fabs(A[j][i]); // ну и по этому циклу 2*(n-i-1) действий
                    i0 = j;
                }
            }
            
            swap(A[i], A[i0]); // O(n) обмен строк, n присваиваний
            swap(b[i], b[i0]); // O(1)

            if (fabs(A[i][i]) < EPS)
            {
                throw std::runtime_error("Система вырождена");
            }


            for (int j = i + 1; j < n; j++)
            {
                real m{ A[j][i] / A[i][i] }; // 1 деление
                for (int k = i; k < n; k++)
                {
                    A[j][k] -= m * A[i][k];  // 1 умн, 1 выч = 2
                } // ну и по этому циклу 2*(n-i) действий
                b[j] -= m * b[i];  // 1 умн, 1 выч = 2
            }  // ну и по этому циклу 2*(n-i-1) + sum(2*(n-i) + 2) = (n-i-1) делений + (n-i-1)*(2*(n-i)+2)
        }


        for (int i = n - 1; i >= 0; i--)
        {
            if (fabs(A[i][i]) < EPS)
            {
                throw std::runtime_error("Система вырождена");
            }
            realS sum{ 0.0 };
            for (int j = i + 1; j < n; j++)
            {
                sum += A[i][j] * x[j]; // 1 умн, 1 сл = 2, тогда 2*(n-i-1)
            }
            x[i] = (b[i] - sum) / A[i][i]; // 1 выч, 1 дел = 2
        } // ну и по этому циклу 2*(n-i-1) + sum(2*(n-i) + 2) = (n-i-1) делений + (n-i-1)*(2*(n-i)+2)
    }
    catch (std::exception& e)
    {
        std::cout << "Ошибка! " << e.what() << std::endl;
        exit(0);
    }
}
void WorkingWithGause()
{
    real** A;
    real* b;
    int n;

    InputMatrix(A, b, n);

    real* x = new real[n];
    for (int i = 0; i < n; i++) x[i] = 0.0;

    SolveGauss(A, b, x, n);

    std::cout << "\nРешение системы методом Гаусса:\n";
    Output(x, n);

    std::cout << "\nМатрица после преобразований:\n";
    PrintMatrix(A, n);

    delete[] b;
    delete[] x;
    for (int i = 0; i < n; i++)
    {
        delete[] A[i];
    }
}
