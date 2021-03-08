#include "SlauSolver.h"

real Matrix::iteration(real* x_k, int i)
{
   real sum = 0;       //
   sum += n - m > i ? Mat[4][i] * x_k[i + m] : 0;  
   sum += n - 1 > i ? Mat[3][i] * x_k[i + 1] : 0;  
   sum +=            Mat[2][i] * x_k[i];                      
   sum += 1 <= i ? Mat[1][i - 1] * x_k[i - 1] : 0; 
   sum += m <= i ? Mat[0][i - m] * x_k[i - m] : 0; 

   return Mat[2][i] ? x_k[i] + w * (b[i] - sum) / Mat[2][i] : 0;
}

void Matrix::Gauss(real **A, real *b, int n, int offset)
{
   Mat = A;
   this->b = b;
   this->n = n;
   m = offset;
   x = new real[n]{};
   int k;
   std::cout.precision(16);
   for (k = 0; k < maxiter; k++)
   {
      for (int i = 0; i < n; i++)
         x[i] = iteration(x, i);

      discr = discrepancy();
      if (discr < eps) break;
      else if (discr >= INFINITY || discr <= -INFINITY)
      {
         std::cout << "\nMatrix is inconsistent\n"; break;
      }
   }
   std::cout << "\nIteration: " << k + 1 << "\tDescrepancy: " << discr << '\n';
   writeSolution();
}

real Matrix::discrepancy()
{
   real sum = 0, sum2 = 0;
   for (int i = 0; i < n; i++)
   {
      sum2 += n - m > i ? Mat[4][i] * x[i + m] : 0;
      sum2 += n - 1 > i ? Mat[3][i] * x[i + 1] : 0;
      sum2 += Mat[2][i] * x[i];
      sum2 += 1 <= i ? Mat[1][i - 1] * x[i - 1] : 0;
      sum2 += m <= i ? Mat[0][i - m] * x[i - m] : 0;

      sum += (b[i] - sum2) * (b[i] - sum2);
      sum2 = 0;
   }
   sum = sqrt(sum);
   return sum / norm(b, n);
}

real Matrix::cond_num(real discr)
{
   real sum = 0;
   real* realx = new real[n];
   for (int i = 0; i < n; i++)
   {
      sum += (i - x[i] + 1) * (i - x[i] + 1);
      realx[i] = i + 1;
   }
   sum = (sqrt(sum) / norm(realx, n)) / discr;
   delete[] realx;
   return sum;
}

real Matrix::norm(real* vec, int n)
{
   real sum = 0; //dub?
   for (int i = 0; i < n; i++)
      sum += (double)vec[i] * (double)vec[i];
   return sqrt(sum);
}

void Matrix::writeSolution()
{
   std::ofstream fS("Solution.txt");
   fS.precision(16);
   for (int i = 0; i < n; i++)
      fS << x[i] << '\n';
}