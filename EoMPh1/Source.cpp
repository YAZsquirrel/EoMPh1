#include "EllepticSolver.h"
#include "SlauSolver.h"

int main()
{
   ElipticEquation Eq;
   Matrix Mat;
   Eq.Init();
   Eq.CreateA();

   for (int t = 0; t < Eq.size; t++) std::cout << Eq.b[t] << '\n';
   std::cout << '\n';

   for (int t = 0; t < Eq.size - Eq.offset;t++) std::cout << Eq.A[0][t] << ' ';
   std::cout << '\n';                 
   for (int t = 0; t < Eq.size - 1;t++) std::cout << Eq.A[1][t]      << ' ';
   std::cout << '\n';                   
   for (int t = 0; t < Eq.size;t++) std::cout << Eq.A[2][t]          << ' ';
   std::cout << '\n';                      
   for (int t = 0; t < Eq.size - 1;t++) std::cout << Eq.A[3][t]          << ' ';
   std::cout << '\n';                     
   for (int t = 0; t < Eq.size - Eq.offset;t++) std::cout << Eq.A[4][t]          << ' ';
   std::cout << '\n';



   Mat.Gauss(Eq.A, Eq.b, Eq.size, Eq.offset);
   Eq.CheckError();

   return 0;
}