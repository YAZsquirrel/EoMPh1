#include "EllepticSolver.h"
#include "SlauSolver.h"

int main()
{
   ElipticEquation Eq;
   Matrix Mat;
   Eq.Init();
   Eq.CreateA();

   //for (int i = 0; i < Eq.size; i++)
   //{
   //    for (int t = 0; t < i - Eq.offset; t++) std::cout << '\t' << 0.;
   //    if (i >= Eq.offset)  std::cout << '\t' << Eq.A[0][i - Eq.offset];
   //
   //    for (int t = 0; i > 1 && t < Eq.offset - 2; t++)
   //        std::cout << '\t' << 0.;
   //
   //    if (i >= 1)	std::cout << '\t' << Eq.A[1][i - 1];
   //    std::cout << '\t' << Eq.A[2][i];
   //     if (i < Eq.size - 1) std::cout << '\t' << Eq.A[3][i]; 
   //
   //    for (int t = 0; t < Eq.offset - 2 && Eq.size - i - 2 > 0; t++) 
   //        std::cout << '\t' << 0.;
   //
   //    if (i < Eq.size - Eq.offset)	std::cout << '\t' << Eq.A[4][i]; 
   //
   //    for (int t = 0; Eq.size - Eq.offset - i - 1> t; t++)
   //        std::cout << '\t' << 0.;
   //    std::cout << '\n';
   //}

   //std::cout << '\n';
   for (int i = 0; i < Eq.size - Eq.offset; i++) std::cout << Eq.A[0][i] << ' '; std::cout << '\n';
   for (int i = 0; i < Eq.size - 1; i++)		 std::cout << Eq.A[1][i] << ' '; std::cout << '\n';
   for (int i = 0; i < Eq.size; i++)			 std::cout << Eq.A[2][i] << ' '; std::cout << '\n';
   for (int i = 0; i < Eq.size - 1; i++)		 std::cout << Eq.A[3][i] << ' '; std::cout << '\n';
   for (int i = 0; i < Eq.size - Eq.offset; i++) std::cout << Eq.A[4][i] << ' '; std::cout << '\n';

   Mat.Gauss(Eq.A, Eq.b, Eq.size, Eq.offset);
   Eq.CheckError();

   return 0;
}