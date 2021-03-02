#include "EllepticSolver.h"

real ElipticEquation::f(int i, int j)
{
   real x1 = x[0] + i * hx;
   real y1 = y[0] + j * hy;
   //return x1 * x1 + y1 * y1 - 4                                                // u^2
   return  x1 * x1 * x1 + y1 * y1 * y1 - 6 * (x1 + y1) ;                         // u^3
   //return x1 * x1 * x1 * x1 + y1 * y1 * y1 * y1 - 12 * (x1 * x1 + y1 * y1);    // u^4
}

void ElipticEquation::UchetKraevyh()
{
   int i, j;
   for (size_t t = 0; t < size; t++)
   {

   i = t % offset; j = t / offset;
   bool left = i == 0;
   bool right = i == nx;
   bool top = j == ny;
   bool bottom = j == 0;

   bool kray1 = left || right || top || bottom;                      //corner
   bool kray2 = !((left || right) && (top || bottom)) && kray1;      //side


   if(kray1 && !kray2)
   {

      if (t < size - offset) A[0][t] =          0.;
      if (t < size - 1)      A[1][t] =          0.;
                             A[2][t] =          1.;
      if (t >= 1)            A[3][t - 1] =      0.;
      if (t >= offset)       A[4][t - offset] = 0.;

      b[t] = ug(i, j);
   }
   else if(kray2){ 
   if (left)
   {

   }
   else if(right)
   {

   }
   else if(top)
   { 
   
   }
   else if (bottom)
   {

   }
   }
   }
}

real ElipticEquation::ug(int i, int j)
{
   real x1 = x[0] + i * hx;
   real y1 = y[0] + j * hy;
   //return x1 * x1 + y1 * y1                               // u^2
   return x1 * x1 * x1 + y1 * y1 * y1;                      // u^3
   //return x1 * x1 * x1 * x1 + y1 * y1 * y1 * y1;          // u^4
}

void ElipticEquation::CheckError()
{
   std::ifstream sol("Solution.txt");
   real sum = 0;
   real buf;
   for (size_t i = 0; i < size; i++)
   {
      sol >> buf;
      buf -= ug(i % offset, i / offset);
      sum += buf * buf;
   }
   std::cout << "\n\nError: " << sqrt(sum) / norm();
}

real ElipticEquation::norm()
{
   real sum = 0;
   for (size_t i = 0; i < size; i++)
      sum += ug(i % offset, i / offset) * ug(i % offset, i / offset);
   return sqrt(sum);
}

void ElipticEquation::CreateA()
{
   int xi, yj;
   for (size_t t = 0; t < size; t++)
   {  
      xi = t % offset; yj = t / offset;

      if (t < size - offset) A[0][t] =          PrimeApproxY();
      if (t < size - 1)      A[1][t] =          PrimeApproxX();
                             A[2][t] =          LaplasApprox();
      if (t >= 1)            A[3][t - 1] =      PrimeApproxX();
      if (t >= offset)       A[4][t - offset] = PrimeApproxY();

      b[t] = f(xi, yj);
   } 

   for (size_t i = 0; i < size; i++)
      std::cout << b[i] << '\n';

#if !DYNAMIC_MESH
   delete[] x;
   delete[] y;
#else // DYNAMIC_MESH
   delete[] hx;
   delete[] hy;
#endif

}

inline real ElipticEquation::PrimeApproxX()
{
   return -lam / (hx * hx);
}

inline real ElipticEquation::PrimeApproxY()
{
   return -lam / (hy * hy);
}

real ElipticEquation::LaplasApprox()
{
   return 2 * lam * (1.f / (hx * hx) + 1.f / (hy * hy)) + gamma;
}

void ElipticEquation::Init()
{
   std::ifstream mesh("Mesh.txt");
   
   mesh >> size >> nx >> ny; 
   offset = nx - 1;
   A = new real*[5];
   b = new real[size]{};
   x = new real[2];
   y = new real[2];
   
   mesh >> x[0] >> x[1];
   mesh >> y[0] >> y[1];

   hx = (x[1] - x[0]) / (nx - 1);
   hy = (y[1] - y[0]) / (ny - 1);

   A[0] = new real[size - offset];
   A[1] = new real[size - 1];
   A[2] = new real[size];
   A[3] = new real[size - 1];
   A[4] = new real[size - offset];
}
