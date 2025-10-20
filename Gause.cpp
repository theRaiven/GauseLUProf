// Gause.cpp

#include "Header.h"
// ����� ������� ������� ������� �� ������
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
                if (fabs(A[j][i]) > maxElem) // 1 ������, 1 ��������� 
                {
                    maxElem = fabs(A[j][i]); // �� � �� ����� ����� 2*(n-i-1) ��������
                    i0 = j;
                }
            }
            
            swap(A[i], A[i0]); // O(n) ����� �����, n ������������
            swap(b[i], b[i0]); // O(1)

            if (fabs(A[i][i]) < EPS)
            {
                throw std::runtime_error("������� ���������");
            }


            for (int j = i + 1; j < n; j++)
            {
                real m{ A[j][i] / A[i][i] }; // 1 �������
                for (int k = i; k < n; k++)
                {
                    A[j][k] -= m * A[i][k];  // 1 ���, 1 ��� = 2
                } // �� � �� ����� ����� 2*(n-i) ��������
                b[j] -= m * b[i];  // 1 ���, 1 ��� = 2
            }  // �� � �� ����� ����� 2*(n-i-1) + sum(2*(n-i) + 2) = (n-i-1) ������� + (n-i-1)*(2*(n-i)+2)
        }


        for (int i = n - 1; i >= 0; i--)
        {
            if (fabs(A[i][i]) < EPS)
            {
                throw std::runtime_error("������� ���������");
            }
            realS sum{ 0.0 };
            for (int j = i + 1; j < n; j++)
            {
                sum += A[i][j] * x[j]; // 1 ���, 1 �� = 2, ����� 2*(n-i-1)
            }
            x[i] = (b[i] - sum) / A[i][i]; // 1 ���, 1 ��� = 2
        } // �� � �� ����� ����� 2*(n-i-1) + sum(2*(n-i) + 2) = (n-i-1) ������� + (n-i-1)*(2*(n-i)+2)
    }
    catch (std::exception& e)
    {
        std::cout << "������! " << e.what() << std::endl;
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

    std::cout << "\n������� ������� ������� ������:\n";
    Output(x, n);

    std::cout << "\n������� ����� ��������������:\n";
    PrintMatrix(A, n);

    delete[] b;
    delete[] x;
    for (int i = 0; i < n; i++)
    {
        delete[] A[i];
    }
}
