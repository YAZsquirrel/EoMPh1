#pragma once
#include <fstream>
#include <vector>
#include <iostream>
typedef double real;

class ElipticEquation
{
private:
   real hx, hy; // шаги

   real x0, x1, y0, y1; //   x0 - левая граница, y0 - нижная, x1 - правая, y1 - верхняя
   real xb, yb;      // границы "вырезанного" прямоугольника
   int nx, ny;		// количество узлов на прямой
   real gamma = 1, lam = 1, theta(int ij, bool xy);
   real q = 2;
   real PrimeApproxXr(int i);	// Приближение производных
   real PrimeApproxYr(int j);
   real PrimeApproxXl(int i);	
   real PrimeApproxYl(int j);
   real LaplasApprox(int i, int j);
   void DivideKnots();
   real f(int i, int j);
   void UchetKraevyh();
   real norm();
public:
   void Init();			// Получаем узлы
   void CreateA();		// собираем матрицу
   real ug(int i, int j);  // по идее - первые краевые, по факту - искомая ф-ия
   void CheckError();      // Смотрим погрешность

   real* b;	// вектор правой части (nx*ny)
   real** A; // глобальная матрица (nx*ny) * (nx*ny) (хотя, учитывая, что всего диагоналей 5, то (5*nx*ny))
   int size = 0;    // колво узлов
   int offset = 0; // количество узлов вдоль оси x
};
