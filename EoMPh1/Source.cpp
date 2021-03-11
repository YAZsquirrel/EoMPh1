#include "EllepticSolver.h"
#include "SlauSolver.h"

int main()
{
   ElipticEquation Eq;
   Matrix Mat;
   Eq.Init();
   Eq.CreateA();

   Mat.Gauss(Eq.A, Eq.b, Eq.size, Eq.offset);
   Eq.CheckError();

   return 0;
}