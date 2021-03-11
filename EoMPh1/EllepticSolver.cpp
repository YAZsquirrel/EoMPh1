#include "EllepticSolver.h"

real ElipticEquation::f(int i, int j)
{
   real x = x0 + i * hx;
   real y = y0 + j * hy;

   //return 1;                                                         // u^0
   //return x + y;                                                     // u^1
   //return x * x + y * y - 4;                                         // u^2
   return x * x * x + y * y * y - 6 * (x + y) ;                      // u^3
   //return x * x * x * x + y * y * y * y - 12 * (x * x + y * y);      // u^4
   return sin(x) + cos(y);
}


real ElipticEquation::ug(int i, int j)
{
   real x = x0 + i * hx;
   real y = y0 + j * hy;

   if (x0 + i * hx > xb && y0 + j * hy > yb)
      return 0.;
   //return 1;                               // u^0
   //return x + y;                           // u^1
   //return x * x + y * y;                   // u^2
   return x * x * x + y * y * y;           // u^3
   //return x * x * x * x + y * y * y * y;   // u^4
   return sin(x) + cos(y);
}

real ElipticEquation::theta(int ij, bool xy)
{
   real th = xy ? x0 + ij * hx : y0 + ij * hy;

   //return 0;                    // u^0
   //return 1; // du/dn           // u^1
   return 2 * th;               // u^2
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
      if (!(i > nxb && j > nyb))
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
               A[2][t] = (lam / hx);
               A[3][t] = (-lam / hx);//-du/dx // t -> t+1 (->)// 
               b[t] = -theta(i, 0);
            }
            else if(right || inner_vert)
            {
               A[2][t] = (lam / hx);
               A[1][t - 1] = (-lam / hx);//du/dx // t -> t-1 (<-) // 

               b[t] = theta(i, 0);
            }
            else if(top || inner_hor)
            { 
               A[2][t] = (lam / hy);
               A[0][t - offset] = (-lam / hy);//du/dy // t -> t - nx (v) // A[0][t] = du/dy (1/hy) 

                b[t] = theta(j, 1);
            }
            else if (bottom)
            {

               A[2][t] = (lam / hy);
               A[4][t] = (-lam / hy);//-du/dy // t -> t + nx (^) // A[4][t] = -du/dy (-1/hy)

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
   std::cout << "\n\nError: " << sqrt(sum) / norm();
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

      outside = x0 + xi * hx > xb && y0 + yj * hy > yb;

      if (t >= offset) A[0][t - offset] = ( outside ? 0. : PrimeApproxY(yj) );
      if (t >= 1)           A[1][t - 1] = ( outside ? 0. : PrimeApproxX(xi) );
                                A[2][t] = ( outside ? 0. : LaplasApprox(xi, yj) );
      if (t < size - 1)         A[3][t] = ( outside ? 0. : PrimeApproxX(xi) );
      if (t < size - offset)    A[4][t] = ( outside ? 0. : PrimeApproxY(yj) );

      b[t] = outside ? 0. : f(xi, yj);
   } 

   UchetKraevyh();

}


inline real ElipticEquation::PrimeApproxX(int i)
{
   return -lam / (hx * hx);
}

inline real ElipticEquation::PrimeApproxY(int j)
{
   return -lam / (hy * hy);
}

real ElipticEquation::LaplasApprox(int i, int j)
{
   return 2 * lam * (1.f / (hx * hx) + 1.f / (hy * hy)) + gamma;
}

void ElipticEquation::DivideKnots()
{
   for (int i = 0; i < size; i++)
   {

   }
   if(q == 1) {
      hx = (x1 - x0) / (nx - 1);
      hy = (y1 - y0) / (ny - 1);
   }
   else {
      hy = (y1 - y0) * (1.f - q) / (1.f - pow(q, nx - 1));
      hx = (x1 - x0) * (1.f - q) / (1.f - pow(q, ny - 1));
   }
}

void ElipticEquation::Init()
{
   std::ifstream mesh("Mesh.txt");
   
   mesh >> nx >> ny;
   size = nx * ny;
   offset = nx;
   A = new real*[5];
   b = new real[size]{};
   
   mesh >> x0 >> x1;
   mesh >> y0 >> y1;
   mesh >> xb >> yb;

   DivideKnots();

   A[0] = new real[size - offset];
   A[1] = new real[size - 1];
   A[2] = new real[size];
   A[3] = new real[size - 1];
   A[4] = new real[size - offset];
}
