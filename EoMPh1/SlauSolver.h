#pragma once
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>

typedef double real;

class Matrix
{
private:
   int n, m, maxiter = 100000;
   real** Mat;
   real* b;
   real eps = 1e-14f, w = 1.51f;
   real* x;
   real discr;
   real discrepancy();       //относительная невязка
   real cond_num(real);			 //число обусловленности
   void writeSolution();
   void read();
   real norm(real* vec, int n);
   real iteration(real* x_k, int i);
public:
   void Gauss(real**, real*, int n, int offset);
};
