// GauseLUProf.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>
#include <fstream>
using namespace std;

using real = float;

void PrintMatrix(real** A, int n)
{
    cout << n << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << A[i][j] << ' ';
        }
        cout << endl;
    }
}
template<class T>
void PrintMatrix(T* A, int n)
{
    for (int j = 0; j < n; j++)
    {
        cout << A[j] << ' ';
    }
    cout << endl;
}

// метод решения в профильном формате
void Solve_Ly_b(const real* di, const real* al, const int* ia, const real* b, real*& y, const int n)
{
    for (int i = 0; i < n; i++)
    {
        real sum{ 0.0 };

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
void Solve_Ux_y(const real* au, const int* ia, const real* y, real*& x, const int n)
{

    for (int i = n - 1; i >= 0; --i)
    {
        real sum{ 0.0 };

        for (int k = i + 1; k < n; ++k)
        {
            int j0{ ia[k] };
            int j1{ ia[k + 1] };
            int j{ k - (j1 - j0) };

            if (i < j) continue;

            int idx{ j0 + (i - j) };
            sum += au[idx] * x[k];
        }

        x[i] = y[i] - sum;
    }
};

void InputProfFormat(real*& di, int*& ia, real*& al, real*& au, real*& b, int& n)
{
    try
    {
        std::ifstream fin("matrixProf.txt");
        if (!fin.is_open())
        {
            throw std::runtime_error("Не удалось открыть файл");
        }

        fin >> n;

        di = new real[n];
        ia = new int[n + 1];
        b = new real[n];

        for (int i = 0; i < n; i++)
        {
            fin >> di[i];
            if (di[i] == 0)
            {
                throw std::runtime_error("Вырожденная система!");
            }
        }
        for (int i = 0; i <= n; i++) fin >> ia[i]; 
        // надо чтобы ia[i] = 0 :)
        int nz = ia[n] - ia[0];

        al = new real[nz];
        au = new real[nz];

        for (int i = 0; i < nz; i++) fin >> al[i];
        for (int i = 0; i < nz; i++) fin >> au[i];

        for (int i = 0; i < n; i++) fin >> b[i];

        fin.close();
    }
    catch (std::exception& e)
    {
        std::cout << "Ошибка! " << e.what() << std::endl;
        exit(0);
    }
};

void ConvertingAtoLUProf(int n, real*& di, int*& ia, real*& al, real*& au)
{
    for (int i = 0; i < n; ++i)
    {
        int i0{ ia[i] };
        int i1{ ia[i + 1] };
        int j_0_i{ i - (i1 - i0) };

        // --- L ---
        for (int k = i0; k < i1; ++k)
        {
            int j{ j_0_i + (k - i0) };

            int j_0_j{ j - (ia[j + 1] - ia[j]) };
            int k0{ max(j_0_i, j_0_j) };
            int k1{ j - 1 };

            real sum{ 0.0 };
            if (k0 <= k1)
            {
                int off_i{ k0 - j_0_i };
                int off_j{ k0 - j_0_j };
                for (int m = 0; m <= k1 - k0; ++m)
                {
                    sum += al[i0 + off_i + m] * au[ia[j] + off_j + m];
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

            real sum{ 0.0 };
            if (k0 <= k1)
            {
                int off_i{ k0 - j_0_i };
                int off_j{ k0 - j_0_j };
                for (int m = 0; m <= k1 - k0; ++m)
                {
                    sum += al[ia[j] + off_j + m] * au[i0 + off_i + m];
                }
            }

            au[k] = (au[k] - sum) / di[j];
        }


        real sumD{ 0.0 };
        for (int k = i0; k < i1; ++k)
        {
            sumD += al[k] * au[k];
        }
        di[i] -= sumD;
    }

    cout << "di: "; PrintMatrix(di, n);
    cout << "ia: "; PrintMatrix(ia, n + 1);
    cout << "al: "; PrintMatrix(al, ia[n]);
    cout << "au: "; PrintMatrix(au, ia[n]);

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
        fout << x[i] << ' ';
        cout << x[i] << ' ';
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
    ConvertingAtoLUProf(n, di, ia, al, au);
    PrintMatrix(b,n);
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

// метод решения плотной матрицы по Гауссу
void InputMatrix(real**& A, real*& b, int& n)
{
    std::ifstream fin("matrixGause.txt");
    if (!fin.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл");
    }

    fin >> n;

    A = new real*[n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new real[n];
    }
    // проверка, что элементы есть
    b = new real[n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fin >> A[i][j];
        }
    }
    for (int i = 0; i < n; i++)
    {
        fin >> b[i];
    }
}

void SolveGauss(real** A, real* b, real*& x, int& n)
{
    try
    {
        const double EPS{ 1e-12 }; // использовать относительный ноль
        for (int i = 0; i < n; i++)
        {
            if (fabs(A[i][i]) < EPS)
            {
                throw std::runtime_error("Система вырождена");
            }


            for (int j = i + 1; j < n; j++)
            {
                double m{ A[j][i] / A[i][i] };

                for (int k = i; k < n; k++)
                {
                    A[j][k] -= m * A[i][k];
                }
                b[j] -= m * b[i];
            }
        }
        // 
        for (int i = n - 1; i >= 0; i--)
        {
            if (fabs(A[i][i]) < EPS)
            {
                throw std::runtime_error("Система вырождена");
            }
            double sum{ 0 };
            for (int j = i + 1; j < n; j++)
            {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        } // swap(a[j], a[i])
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
    real* x = new real[n] ;
    SolveGauss(A, b, x, n);
    Output(x, n);

    delete[] b;
    delete[] x;
    for (int i = 0; i < n; i++)
    {
        delete[] A[i];
    }
}

// Конвертор

void ConvertingMatrixAtoLU(real** A, int n, real*& di, int*& ia, real*& al, real*& au)
{
    real** L = new real * [n];
    real** U = new real * [n];
    for (int i = 0; i < n; i++)
    {
        L[i] = new real[n]{ 0 };
        U[i] = new real[n]{ 0 };
    }
    for (int i = 0; i < n; i++)
    {
        U[i][i] = 1;
        L[i][i] = A[i][i];

        for (int j = 0; j <= i; j++)
        {
            float sum{ 0 };
            for (int k = 0; k < j; k++)
            {
                sum += L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - sum;

        }
        for (int j = n - 1; j > i; j--)
        {
            float sum{ 0 };
            for (int k = 0; k < i; k++)
            {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = (A[i][j] - sum) / L[i][i];
        }
    }
    cout << "L:\n"; PrintMatrix(L, n); cout << "\nU:\n"; PrintMatrix(U, n);
    di = new real[n];
    ia = new int[n + 1];
    int nz = 0;
    ia[0] = 0;
    
    // Освобождаем временные матрицы
    for (int i = 0; i < n; i++)
    {
        delete[] L[i];
        delete[] U[i];
    }
    delete[] L;
    delete[] U;
}
void WorkingWithConverter()
{
    real** A;
    real* di; int* ia;
    real* al; real* au;
    real* b;
    int n;
    InputMatrix(A, b, n);
    ConvertingMatrixAtoLU(A, n, di, ia, al, au);
    //InputProfFormat(di, ia, al, au, b, n);

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
    for (int i = 0; i < n; i++)
    {
        delete[] A[i];
    }
}


int main()
{
    setlocale(LC_ALL, "rus");

    //WorkingWithGause();
    WorkingWithProfFormat();
    //WorkingWithConverter();

}
