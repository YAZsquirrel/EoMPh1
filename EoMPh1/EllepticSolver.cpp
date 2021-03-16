#include "EllepticSolver.h"

real ElipticEquation::f(int i, int j)
{
   real x = x0;
   real y = y0;
   for (size_t p = 0; p < i; p++)
       x += hx * pow(q, p);
   for (size_t p = 0; p < j; p++)
       y += hy * pow(q, p);
   if (x > xb && y > yb)
      return 0.;
   //return 1;                                                         // u^0
   //return x + y;                                                     // u^1
   return x * x + y * y - 4;                                         // u^2
   //return x * x * x + y * y * y - 6 * (x + y) ;                      // u^3
   //return x * x * x * x + y * y * y * y - 12 * (x * x + y * y);      // u^4
   return sin(x) + cos(y);
}

real ElipticEquation::ug(int i, int j)
{
   real x = x0;
   real y = y0;
   for (size_t p = 0; p < i; p++)
      x += hx * pow(q, p);
   for (size_t p = 0; p < j; p++)
      y += hy * pow(q, p);

   if (x > xb && y > yb)
      return 0.;
   //return 1;                               // u^0
   //return x + y;                           // u^1
   return x * x + y * y;                   // u^2
   //return x * x * x + y * y * y;           // u^3
   //return x * x * x * x + y * y * y * y;   // u^4
   return sin(x) + cos(y);
}

real ElipticEquation::theta(int ij, bool xy)
{
   real th = xy ? x0 + ij * hx * pow(q, ij) : y0 + ij * hy * pow(q, ij);

   //return 0;                    // u^0
   return 1; // du/dn           // u^1
   //return 2 * th;               // u^2
   //return 3 * th * th;          // u^3
   //return 4 * th * th * th;     // u^4
}

void ElipticEquation::UchetKraevyh()
{
   int i, j;
   int nyb = (int)(ny * yb / (y1 - y0));
   int nxb = (int)(nx * xb / (x1 - x0));
   bool left, right, top, bottom;
   bool inner_vert, inner_hor; // |_
   bool kray1, kray2;
   
   for (int t = 0; t < size; t++)
   {
      i = t % nx; j = t / nx;

      real x = x0;
      real y = y0;
      for (size_t p = 0; p < i; p++)
         x += hx * pow(q, p);
      for (size_t p = 0; p < j; p++)
         y += hy * pow(q, p);

      if (!(x > xb && y > yb))
      {
         left = i == 0; right = i == nx - 1;
         top = j == ny - 1; bottom = j == 0;
         inner_hor = j == nyb && i >= nxb; inner_vert = i == nxb && j >= nyb;

         kray1 = left || right || top || bottom || ((inner_hor || inner_vert) && !(inner_hor && inner_vert)) ;           //corner
         kray2 = !(((left || right) && (top || bottom)) || ((right && inner_hor) || (top && inner_vert))) && kray1;      //side
            //  не на углах большого прямоугольника(без выреза) и не на углах, вырезанные малым прямоугольником, но на краях 

         if(kray1)
         {     
            if (t >= offset) A[0][t - offset] = 0.;
            if (t >= 1)           A[1][t - 1] = 0.;
                                      A[2][t] = 1.;
            if (t < size - 1)         A[3][t] = 0.;
            if (t  < size - offset)   A[4][t] = 0.;
         
            b[t] = ug(i, j);
         }

         if(kray2 && false){
            if (left)
            {
               A[2][t] = (lam / hx * pow(q, i));
               A[3][t] = (-lam / hx * pow(q, i));//-du/dx // t -> t+1 (->)// 
               b[t] = -theta(i, 0);
            }
            else if(right || inner_vert)
            {
               A[2][t] = (lam / hx * pow(q, i - 1));
               A[1][t - 1] = (-lam / hx * pow(q, i - 1));//du/dx // t -> t-1 (<-) // 

               b[t] = theta(i, 0);
            }
            else if(top || inner_hor)
            { 
               A[2][t] = (lam / hy * pow(q, j - 1));
               A[0][t - offset] = (-lam / hy * pow(q, j - 1));//du/dy // t -> t - nx (v) // A[0][t] = du/dy (1/hy) 

                b[t] = theta(j, 1);
            }
            else if (bottom)
            {

               A[2][t] = (lam / hy * pow(q, j));
               A[4][t] = (-lam / hy * pow(q, j));//-du/dy // t -> t + nx (^) // A[4][t] = -du/dy (-1/hy)

                b[t] = -theta(j, 1);
            }
         }         
      }
   }
}

void ElipticEquation::CheckError()
{
   std::ifstream sol("Solution.txt");
   real sum = 0;
   real buf = 0;
   for (int i = 0; i < size; i++)
   {
      sol >> buf;
      std::cout << buf << " - " << ug(i % nx, i / nx) << '\n';
      buf -= ug(i % nx, i / nx);
      sum += buf * buf;
   }
   std::cout.precision(7);
   std::cout << "\n\nError: " << std::scientific << sqrt(sum) / norm();
}

real ElipticEquation::norm()
{
   real sum = 0;
   for (int i = 0; i < size; i++)
      sum += ug(i % offset, i / offset) * ug(i % offset, i / offset);
   return sqrt(sum);
}

void ElipticEquation::CreateA()
{
   int xi, yj;
   bool outside;
   for (int t = 0; t < size; t++)
   {  
      xi = t % nx; yj = t / nx;

      real x = x0;
      real y = y0;
      for (int p = 0; p < xi; p++)
         x += hx * pow(q, p);
      for (int p = 0; p < yj; p++)
         y += hy * pow(q, p);

      outside = x > xb && y > yb; //xi yj?
      outside = outside;
      if (t >= offset) A[0][t - offset] = ( outside ? 0. : PrimeApproxYl(yj) );
      if (t >= 1)           A[1][t - 1] = ( outside ? 0. : PrimeApproxXl(xi) );
                                A[2][t] = ( outside ? 0. : LaplasApprox(xi, yj) );
      if (t < size - 1)         A[3][t] = ( outside ? 0. : PrimeApproxXr(xi) );
      if (t < size - offset)    A[4][t] = ( outside ? 0. : PrimeApproxYr(yj) );

      b[t] = outside ? 0. : f(xi, yj);
   }
   UchetKraevyh();

}

inline real ElipticEquation::PrimeApproxXr(int i)
{
   return -2.f * lam / (hx * pow(q, i) * hx * (pow(q, i) + pow(q, i-1)));
}

inline real ElipticEquation::PrimeApproxYr(int j)
{
   return -2.f * lam / (hy * pow(q, j) * hy * (pow(q, j) + pow(q, j-1)));
}

inline real ElipticEquation::PrimeApproxXl(int i)
{
   return -2.f * lam / (hx * pow(q, i-1) * hx * (pow(q, i) + pow(q, i-1)));
}

inline real ElipticEquation::PrimeApproxYl(int j)
{
   return -2.f * lam / (hy * pow(q, j-1) * hy * (pow(q, j) + pow(q, j-1)));
}

real ElipticEquation::LaplasApprox(int i, int j)
{
   return 2.f * lam * (1.f / (hx * pow(q, 2 * i - 1) * hx) + 1.f / (hy * pow(q, 2 * j - 1) * hy)) + gamma;
}

void ElipticEquation::DivideKnots()
{
   if(q == 1.f) {
      hx = (x1 - x0) / (nx - 1);
      hy = (y1 - y0) / (ny - 1);
   }
   else {
      hx = (x1 - x0) * (1.f - q) / (1.f - pow(q, nx - 1));
      hy = (y1 - y0) * (1.f - q) / (1.f - pow(q, ny - 1));
      real sum = x0;
      for (int i = 0; i < nx - 1; i++)
      {
         std::cout << sum << " + " << hx * pow(q, i) << '\n';
         sum += hx * pow(q, i);
      }
      std::cout << sum << " + 0\n";
   }                                                         
}

void ElipticEquation::Init()
{
   std::ifstream mesh("Mesh.txt");
   
   mesh >> q >> nx >> ny;
   size = nx * ny;
   offset = nx;
   A = new real*[5];
   b = new real[size]{};
   
   mesh >> x0 >> x1;
   mesh >> y0 >> y1;
   mesh >> xb >> yb;

   DivideKnots();

   A[0] = new real[size - offset]{};
   A[1] = new real[size - 1]{};
   A[2] = new real[size]{};
   A[3] = new real[size - 1]{};
   A[4] = new real[size - offset]{};
}
