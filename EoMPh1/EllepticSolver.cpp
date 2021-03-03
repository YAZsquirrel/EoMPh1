#include "EllepticSolver.h"

real ElipticEquation::f(int i, int j)
{
   real x1 = x0 + i * hx;
   real y1 = y0 + j * hy;
   //return x1 * x1 + y1 * y1 - 4;                                                // u^2
   //return  x1 * x1 * x1 + y1 * y1 * y1 - 6 * (x1 + y1) ;                       // u^3
   return x1 * x1 * x1 * x1 + y1 * y1 * y1 * y1 - 12 * (x1 * x1 + y1 * y1);    // u^4
}

void ElipticEquation::UchetKraevyh()
{
   int i, j;
   int nyb = (int)(ny * yb / (y1 - y0));
   int nxb = (int)(nx * xb / (x1 - x0));
   bool left, right, top, bottom;
   bool inner_vert, inner_hor; // |_
   bool kray1, kray2;
   
   for (size_t t = 0; t < size; t++)
   {
      i = t % nx; j = t / nx;
      if (!(i > nxb && j > nyb))
      {
         left = i == 0; right = i == nx - 1;
         top = j == ny - 1; bottom = j == 0;
         inner_hor = j == nyb && i >= nxb; inner_vert = i == nxb && j >= nyb;

         kray1 = left || right || top || bottom || ((inner_hor || inner_vert) && !(inner_hor && inner_vert)) ;                      //corner
         kray2 = !(((left || right) && (top || bottom)) || ((right && inner_hor) || (top && inner_vert))) && kray1;      //side
            //  не на углах большого прямоугольника(без выреза) и не на углах, вырезанные малым треугольником, но на краях 

         if(kray1)
         {

            if (t < size - offset) A[0][t] =          0.;
            if (t < size - 1)      A[1][t] =          0.;
                                   A[2][t] =          1.;
            if (t >= 1)            A[3][t - 1] =      0.;
            if (t >= offset)       A[4][t - offset] = 0.;

            b[t] = ug(i, j);
         }
         if(kray2){ 
            if (left)
            {
                A[3][t + 1] = (-1 / hx);//-du/dx // t -> t+1 (->)// 
            }
            else if(right || inner_vert)
            {
                A[1][t - 1] = (1 / hx);//du/dx // t -> t-1 (<-) // 
            }
            else if(top || inner_hor)
            { 
                A[0][t] = (1 / hy);//du/dy // t -> t - nx (v)// A[0][t] = du/dy (1/hy) 
            }
            else if (bottom)
            {
                A[4][t] = (-1 / hy);//-du/dy // t -> t + nx (^)// A[4][t] = -du/dy (-1/hy) 
            }
         }
      }
   }
}

real ElipticEquation::ug(int i, int j)
{
   if (x0 + i * hx > xb && y0 + j * hy > yb)
      return 0.;
   real x1 = x0 + i * hx;
   real y1 = y0 + j * hy;
   //return x1 * x1 + y1 * y1;                              // u^2
   //return x1 * x1 * x1 + y1 * y1 * y1;                    // u^3
   return x1 * x1 * x1 * x1 + y1 * y1 * y1 * y1;          // u^4
}

void ElipticEquation::CheckError()
{
   std::ifstream sol("Solution.txt");
   real sum = 0;
   real buf = 0;
   for (size_t i = 0; i < size; i++)
   {
      sol >> buf;
      buf -= ug(i % nx, i / nx);
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
   bool outside;
   for (size_t t = 0; t < size; t++)
   {  
      xi = t % nx; yj = t / nx;

      outside = x0 + xi * hx > xb && y0 + yj * hy > yb;

      if (t < size - offset) A[0][t] =          ( outside ? 0. : PrimeApproxY() );
      if (t < size - 1)      A[1][t] =          ( outside ? 0. : PrimeApproxX() );
                             A[2][t] =          ( outside ? 0. : LaplasApprox() );
      if (t >= 1)            A[3][t - 1] =      ( outside ? 0. : PrimeApproxX() );
      if (t >= offset)       A[4][t - offset] = ( outside ? 0. : PrimeApproxY() );

      b[t] = outside ? 0. : f(xi, yj);
   } 

   UchetKraevyh();

   for (size_t i = 0; i < size; i++)
      std::cout << b[i] << '\n';

#if !DYNAMIC_MESH
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
   offset = nx;
   A = new real*[5];
   b = new real[size]{};
   
   mesh >> x0 >> x1;
   mesh >> y0 >> y1;
   mesh >> xb >> yb;

   hx = (x1 - x0) / (nx - 1);
   hy = (y1 - y0) / (ny - 1);

   A[0] = new real[size - offset];
   A[1] = new real[size - 1];
   A[2] = new real[size];
   A[3] = new real[size - 1];
   A[4] = new real[size - offset];
}
