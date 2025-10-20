// Header.h

#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

using real = double;
using realS = double;
const int accuracy = 15;

//
template<typename T>
constexpr T relativeEPS();

template<>
constexpr float relativeEPS<float>() { return 1e-34f; }

template<>
constexpr double relativeEPS<double>() { return 1e-300; }
//
template<class T>
void PrintMatrix(T** A, int n)
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
void PrintMatrix(T* A, int n);

// метод решения в профильном формате
void Solve_Ly_b(const real* di, const real* al, const int* ia, const real* b, real*& y, const int n);
void Solve_Ux_y(const real* au, const int* ia, real* y, real*& x, const int n);
template<class T>
void InputProfFormat(real*& di, int*& ia, real*& al, real*& au, T*& b, int& n)
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
        b = new T[n];

        for (int i = 0; i < n; i++)
        {
            fin >> di[i];
            if (di[i] == 0)
            {
                throw std::runtime_error("Вырожденная система!");
            }
        }
        for (int i = 0; i <= n; i++) fin >> ia[i];
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

void ConvertingAtoLUProf( real*& di, int*& ia, real*& al, real*& au, int n);
void Output(real* x, int n);
void WorkingWithProfFormat();
// метод решения плотной матрицы по Гауссу
template<class T>
void InputMatrix(real**& A, T*& b, int& n)
{
    std::ifstream fin("matrixGause.txt");
    if (!fin.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл");
    }

    fin >> n;

    A = new real * [n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new real[n];
    }
    b = new T[n];
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


void SolveGauss(real** A, real* b, real*& x, int& n);
void WorkingWithGause();

// Решение слау A^k * x^k = F^k
void MultDense(const real** A, const int* xTrue, real* b, int n);
void WorlingWithDense();
void SubtractVectors(int n, const int* xTrue, const real* xFound);


void WorkingWithHilbert();
void MultSkyline(const real* di, const int* ia, const real* al, const real* au, const int* xTrue, real* b, int n);