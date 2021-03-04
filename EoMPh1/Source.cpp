#include "EllepticSolver.h"
#include "SlauSolver.h"

int main()
{
   ElipticEquation Eq;
   Matrix Mat;
   Eq.Init();
   Eq.CreateA();

   std::cout << '\n';
   for (int i = 0; i < Eq.size - Eq.offset; i++) std::cout << Eq.A[0][i] << ' '; std::cout << '\n';
   for (int i = 0; i < Eq.size - 1; i++)		 std::cout << Eq.A[1][i] << ' '; std::cout << '\n';
   for (int i = 0; i < Eq.size; i++)			 std::cout << Eq.A[2][i] << ' '; std::cout << '\n';
   for (int i = 0; i < Eq.size - 1; i++)		 std::cout << Eq.A[3][i] << ' '; std::cout << '\n';
   for (int i = 0; i < Eq.size - Eq.offset; i++) std::cout << Eq.A[4][i] << ' '; std::cout << '\n';

   Mat.Gauss(Eq.A, Eq.b, Eq.size, Eq.offset);
   Eq.CheckError();

   return 0;
}