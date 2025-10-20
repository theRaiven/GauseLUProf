// LU.cpp

#include "Header.h"

template<class T>
void PrintMatrix(T* A, int n)
{
    for (int j = 0; j < n; j++)
    {
        cout << std::setprecision(accuracy) << A[j] << ' ';
    }
    cout << endl;

}


// метод решения в профильном формате
void Solve_Ly_b(const real* di, const real* al, const int* ia, const real* b, real*& y, const int n)
{
    for (int i = 0; i < n; i++)
    {
        realS sum{ 0.0 };

        int i0{ ia[i] };
        int i1{ ia[i + 1] };

        int j{ i - (i1 - i0) };

        for (int k = i0; k < i1; k++, j++)
        {
            sum += al[k] * y[j];
        }

        y[i] = (b[i] - sum) / di[i];
    }
};
void Solve_Ux_y(const real* au, const int* ia, real* y, real*& x, const int n)
{

    for (int i = n - 1; i >= 0; --i)
    {
        //realS sum{ 0.0 };

        for (int k = i + 1; k < n; ++k)
        {
            int j0{ ia[k] };
            int j1{ ia[k + 1] };
            int j{ k - (j1 - j0) };

            if (i < j) continue;

            int idx{ j0 + (i - j) };
            y[i] -= au[idx] * y[k];
        }
        x[i] = y[i];
    }
};


void ConvertingAtoLUProf(real*& di, int*& ia, real*& al, real*& au, int n)
{
    for (int i = 0; i < n; i++)
    {
        int i0{ ia[i] };
        int i1{ ia[i + 1] };
        int j_0_i{ i - (i1 - i0) };

        // --- L ---
        for (int k = i0; k < i1; k++)
        {
            int j{ j_0_i + (k - i0) };

            int j_0_j{ j - (ia[j + 1] - ia[j]) };
            int k0{ max(j_0_i, j_0_j) };
            int k1{ j - 1 };

            realS sum{ 0.0 };
            if (k0 <= k1)
            {
                int off_i{ k0 - j_0_i };
                int off_j{ k0 - j_0_j };
                for (int m = 0; m <= k1 - k0; m++)
                {
                    sum += (al[i0 + off_i + m]) * (au[ia[j] + off_j + m]);
                }
            }

            al[k] -= sum;
        }

        // --- U ---
        for (int k = i0; k < i1; ++k)
        {
            int j{ j_0_i + (k - i0) };

            int j_0_j{ j - (ia[j + 1] - ia[j]) };
            int k0{ max(j_0_i, j_0_j) };
            int k1{ j - 1 };

            realS sumS{ 0.0 };
            if (k0 <= k1)
            {
                int off_i{ k0 - j_0_i };
                int off_j{ k0 - j_0_j };
                for (int m = 0; m <= k1 - k0; ++m)
                {
                    sumS += (al[ia[j] + off_j + m]) * (au[i0 + off_i + m]);

                }
            }

            au[k] = (au[k] - sumS) / di[j];
        }


        realS sumD{ 0.0 };
        for (int k = i0; k < i1; ++k)
        {
            sumD += (al[k] * au[k]);
        }
        di[i] -= sumD;
    }

    /*cout << "di: "; PrintMatrix(di, n);
    cout << "ia: "; PrintMatrix(ia, n + 1);
    cout << "al: "; PrintMatrix(al, ia[n]);
    cout << "au: "; PrintMatrix(au, ia[n]);*/

}
void Output(real* x, int n)
{
    std::ofstream fout("vectorX.txt");

    if (!fout.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл");
    }
    for (int i = 0; i < n; i++)
    {
        fout << std::setprecision(accuracy) << x[i] << ' ';
        cout << std::setprecision(accuracy) << x[i] << endl;
    }
    cout << endl;
    fout.close();
}
void WorkingWithProfFormat()
{
    int n;

    real* di; int* ia;
    real* al; real* au;
    real* b;

    InputProfFormat(di, ia, al, au, b, n);
    ConvertingAtoLUProf( di, ia, al, au, n);
    PrintMatrix(b, n);
    real* y = new real[n];
    real* x = new real[n];

    Solve_Ly_b(di, al, ia, b, y, n);
    Solve_Ux_y(au, ia, y, x, n);

    Output(x, n);

    delete[] di;
    delete[] ia;
    delete[] al;
    delete[] au;
    delete[] b;
    delete[] x;
    delete[] y;
}