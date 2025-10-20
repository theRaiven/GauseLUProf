// HilbertMult.cpp

#include "Header.h"

void HilbertMatrix(int n, double**& A)
{
    A = new double*[n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            A[i][j] = 1.0 / (i + j + 1);
        }
    }
}

void ConvertToProfile(double** A, int n,
    real*& di, int*& ia, real*& al, real*& au)
{
    di = new real[n];
    for (int i = 0; i < n; i++) di[i] = A[i][i];

    ia = new int[n + 1];
    ia[0] = 0;
    int nz{ 0 };

    for (int i = 0; i < n; i++)
    {
        int j0{ i };
        for (int j = 0; j < i; j++)
        {
            if (A[i][j] != 0.0)
            {
                j0 = j; 
                break;
            }
        }
        int len{ i - j0 };
        ia[i + 1] = ia[i] + len;
        nz += len;
    }

    al = new real[nz];
    au = new real[nz];

    for (int i = 0; i < n; i++)
    {
        int i0{ ia[i] };
        int i1{ ia[i + 1] };
        int len { i1 - i0 };
        for (int k = i0; k < i1; k++) 
        {
            int j = i - i1 + k;
            al[k] = A[i][j];
            au[k] = A[j][i];
        }
    }
}
void MultSkyline(const real* di, const int* ia, const real* al, const real* au, const int* xTrue, double* b, int n)
{
    for (int i = 0; i < n; i++) b[i] = 0.0;
    for (int i = 0; i < n; i++)
    {
        b[i] += di[i] * xTrue[i];
    }
    for (int i = 1; i < n; i++)
    {
        int i0 = ia[i];
        int i1 = ia[i + 1];
        int j0 = i - (i1 - i0);
        int k = i0;
        for (int j = j0; j < i; j++, k++)
        {
            b[i] += al[k] * xTrue[j];
            b[j] += au[k] * xTrue[i];
        }
    }
}
void WorkingWithHilbert()
{
    string filename = "matrixProf.txt";

    for (int n = 1; n <= 14; n++)
    {
        cout << "\n=== Матрица Гильберта " << n << "x" << n << " ===" << endl;

        double** A;

        HilbertMatrix(n, A);
        PrintMatrix(A, n);
        real* di; int* ia; real* al; real* au;
        ConvertToProfile(A, n, di, ia, al, au);

        int* xTrue = new int[n];
        for (int i = 0; i < n; i++) xTrue[i] = i + 1;

        double* b = new double[n];
        MultSkyline(di, ia, al, au, xTrue, b, n);

        ConvertingAtoLUProf(di, ia, al, au, n);

        real* y = new real[n];
        real* x = new real[n];
        real* bReal = new real[n];
        for (int i = 0; i < n; i++)
        {
            bReal[i] = static_cast<real>(b[i]);
        }

        Solve_Ly_b(di, al, ia, bReal, y, n);
        Solve_Ux_y(au, ia, y, x, n);
        SubtractVectors(n, xTrue, x);

        cout << "\nx^k:" << endl;
        Output(x, n);

        delete[] di; delete[] ia; delete[] al; delete[] au; delete[] xTrue; delete[] b;
        delete[] y; delete[] x;
    }
}