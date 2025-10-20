// main.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <format>
#include "Header.h"

int main()
{
    setlocale(LC_ALL, "rus");

    float a = 15;
    int k = 5;
    float a11 = a + pow(10, -k);
    cout << format("{:.7e}",a11);
    //WorlingWithDense();
}

