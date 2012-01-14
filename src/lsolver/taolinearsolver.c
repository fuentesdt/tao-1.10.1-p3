#include "tao_solver.h"    /*I "tao_solver.h"  I*/
#include "taolinearsolver.h"

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolver::PreSolve"
/*@C
   PreSolve - Prepares a matrix for the linear solve.

   Input Parameter:
.  M - the matrix

   Level: intermediate
@*/
int TaoLinearSolver::PreSolve(TaoMat *M)
{
  TaoFunctionBegin;
  SETERRQ(56, "Operation not defined");
  TaoFunctionReturn(1);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolver::Solve"
/*@C
   Solve - Finds the solution to A x = b .

   Input Parameter:
.  b - the right-hand side

   Output Parameter:
+  x - the solution
-  flag - TAO_TRUE if successful solve, TAO_FALSE otherwise.

   Level: intermediate
@*/
int TaoLinearSolver::Solve(TaoVec *b, TaoVec *x, TaoTruth *flag)
{
  TaoFunctionBegin;
  SETERRQ(56, "Operation not defined");
  TaoFunctionReturn(1);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::SolveTrustRegion"
/*@C
   SolveTrustRegion - Finds the solution to the trust-region problem 
       min <x, Ax - b>
       st  |x| <= r
   If A is a symmetrix positive definite matrix and r is large enough,
   then x solves the system of equations Ax = b.
   
   Input Parameter:
.  b - the right-hand side
.  r - the trust-region radius

   Output Parameter:
+  x - the solution
-  flag - TAO_TRUE if successful solve, TAO_FALSE otherwise.

   Level: intermediate
@*/
int TaoLinearSolver::SolveTrustRegion(TaoVec *b, TaoVec *x, double r, 
                                      TaoTruth *t)
{
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  TaoFunctionReturn(1);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolver::GetNumberIterations"
/*@C
   GetNumberIterations - Gets the number of iterations used to solve the system.

   Output Parameter:
.  iters - the number of iterations

   Level: intermediate
@*/
int TaoLinearSolver::GetNumberIterations(int * iters){
  int info;
  TaoFunctionBegin;
  info = PetscInfo(0,"TaoLinearSolver::GetNumberIterations() iterations not available.\n"); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolver::SetTolerances"
/*@C 
  SetTolerances - Sets the convergence tolerances for solving the linear system
  
  Input Parameters:
+ rtol - the relative tolerance
. atol - the absolute tolerance
. dtol - the divergence tolerance
- maxits - the maximum number of iterates

  Level: intermediate
@*/  

int TaoLinearSolver::SetTolerances(double rtol, double atol, double dtol, int maxits) {
  int info;
  info = PetscInfo(0,"inherited method not defined.  Tolerances ignored.\n");
  CHKERRQ(info);
  TaoFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolver::SetOptions"
/*@C
  SetOptions - Sets any relevant options from the command line.

  Input Parameters: none
  
  Level: intermediate
@*/
int TaoLinearSolver::SetOptions(){
  int info;
  TaoFunctionBegin;
  info = PetscInfo(0,"TaoLinearSolver::SetOptions(): method not defined or implemented.\n"); CHKERRQ(info);
  info = TaoOptionsHead("TaoLinearSolverOptions:  Not currently implemented through TAO");CHKERRQ(info);
  info = TaoOptionsTail();CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolver::View"
/*@C
  View - Views the linear solver.

  Input Parameters: none
  
  Level: intermediate
@*/
int TaoLinearSolver::View(){
  int info;
  TaoFunctionBegin;
  info = PetscInfo(0,"TaoLinearSolver::View() method not defined.\n"); CHKERRQ(info);
  info = TaoPrintStatement(0,"  TaoLinearSolverView: No information available: \n");CHKERRQ(info);
  TaoFunctionReturn(0);
}

