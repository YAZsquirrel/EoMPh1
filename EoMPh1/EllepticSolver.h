#pragma once
#include <fstream>
#include <vector>
#include <iostream>
#define DYNAMIC_MESH false
typedef double real;

class ElipticEquation
{
private:
#if !DYNAMIC_MESH
   real hx, hy; // шаги; статичная сетка
#else // DYNAMIC_MESH
   real *hx, *hy; // шаги; динамическая сетка (nx, ny)
#endif
   real *x, *y; //???
   real **xy;
   size_t nx, ny;		// количество узлов на прямой
   std::vector<real> mesh;
   real gamma = 1, lam = 1;
   real PrimeApproxX();	// Приближение производных
   real PrimeApproxY();
   real LaplasApprox();
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
   real** A; // глобальная матрица (nx*ny)*(nx*ny) (хотя, учитывая, что всего диагоналей 5, то (5*nx*ny))
   size_t size = 0;    // колво узлов
   size_t offset = 0; // количество узлов вдоль оси x
};
