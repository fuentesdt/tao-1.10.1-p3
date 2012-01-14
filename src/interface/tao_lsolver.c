#include "src/tao_impl.h"      /*I "tao_solver.h"  I*/


#undef __FUNCT__  
#define __FUNCT__ "TaoCreateLinearSolver"
/*@C
  TaoCreateLinearSolver - Create a linear solver for use in TAO

   Input Parameters:
+  tao - the TAO_SOLVER context
.  MM  - the matrix associated with the solver
-  stype - the type of linear solver

   Output Parameters:
.  ksp - a linear solver

   Types of linear solvers:
+  100 - Solver for any square matrix
.  110 - GMRES
.  200 - Any symmetric matrix
.  210 - MINRES
.  220 - CG with Trust Region
.  230 - SYMMLQ
.  300 - Any symmetric positive definite
-  310 - Conjugate Gradient

   Level: developer

.seealso: TaoLinearSolve()

.keywords: linear solver
@*/
int TaoCreateLinearSolver(TAO_SOLVER tao, TaoMat *MM, TaoInt stype, TaoLinearSolver **ksp)
{
  int info;
  TaoVec *xx;
  TaoLinearSolver* nksp=0;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = TaoGetSolution(tao,&xx);CHKERRQ(info);
  if (MM){
    if (tao->taoappl){
      info = tao->taoappl->GetLinearSolver(MM,stype,&nksp);CHKERRQ(info);
      if (ksp){
	*ksp = nksp;
      } 
      tao->ksp = nksp;
    }
  } else {
    SETERRQ(1,"No matrix specified.");
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoDestroyLinearSolver"
/*@C
  TaoDestroyLinearSolver - Destroy the linear solver used in TAO

   Input Parameters:
.  tao - the TAO_SOLVER context

   Level: developer

.seealso: TaoGetLinearSolver()

.keywords: linear solver
@*/
int TaoDestroyLinearSolver(TAO_SOLVER tao)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  tao->ksp=0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoLinearSolve"
/*@C
  TaoLinearSolve - Solve a linear system

   Input Parameters:
+  tao - the TAO_SOLVER context
.  AA - the matrix
-  bb - the right hand side

   Output Parameters:
+  xx - the solution
-  success - false if linear solve is not successful 

   Level: developer

.seealso: TaoSolve()

.keywords: linear solver
@*/
int TaoLinearSolve(TAO_SOLVER tao, TaoMat *AA, TaoVec *bb, TaoVec* xx, TaoTruth *success)
{
  int info,lits;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (tao->ksp){
    info = tao->ksp->Solve(bb,xx,success);CHKERRQ(info);
    info = tao->ksp->GetNumberIterations(&lits);CHKERRQ(info);
    tao->linear_its += lits;
    TaoFunctionReturn(0);
  } else {
    TaoFunctionReturn(1);
  }
}

#undef __FUNCT__  
#define __FUNCT__ "TaoLinearSolveTrustRegion"
/*@C
  TaoLinearSolveTrustRegion - Minimize the quadratic function <x, Qx - b> 
  with the trust region constraint norm(x, M) <= r, where M is a symmetric 
  positive definite preconditioner.

   Input Parameters:
+  tao - the TAO_SOLVER context
.  Q - the matrix
.  b - the right hand side
-  r - the trust region radius

   Output Parameters:
+  x - the solution
-  success - false if minimization is not successful

   Level: developer

.seealso: TaoLinearSolve()

.keywords: linear solver
@*/
int TaoLinearSolveTrustRegion(TAO_SOLVER tao, TaoMat *Q, TaoVec *b, TaoVec *x, 
			      double r, TaoTruth *success)
{
  int lits,info;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = tao->ksp->SolveTrustRegion(b, x, r, success); CHKERRQ(info);
  info = tao->ksp->GetNumberIterations(&lits);CHKERRQ(info);
  tao->linear_its += lits;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPreLinearSolve"
/*@C
  TaoPreLinearSolve - Prepare to solve a linear system

   Input Parameters:
+  tao - the TAO_SOLVER context
-  AA - the matrix

   Level: developer

.seealso: TaoLinearSolve()

.keywords: linear solver
@*/
int TaoPreLinearSolve(TAO_SOLVER tao, TaoMat *AA)
{
  int info;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = tao->ksp->PreSolve(AA);CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetLinearSolver"
/*@C
   TaoGetLinearSolver - Returns the KSP context for a TAO_SOLVER solver.

   Not Collective, but if TAO_SOLVER object is parallel, then KSP object is parallel

   Input Parameter:
.  solver -  a TAO optimization solver

   Output Parameter:
.  ksp - the KSP context

   Notes:
   The user can then directly modify the linear solver context to modify the Krylov method, preconditioner,
   and tolerances.

   Level: developer

.keywords: Linear Solver, context
x
.seealso: TaoGetKSP()
@*/
int TaoGetLinearSolver(TAO_SOLVER tao,TaoLinearSolver **ksp)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (ksp){
    *ksp = tao->ksp;
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoViewLinearSolver"
/*@C
  TaoViewLinearSolver - View the linear solver

   Input Parameters:
.  tao - the TAO_SOLVER context

   Options Database Key:
.  -tao_kspview - Calls TaoViewLinearSolver() at end of TaoSolve()

   Level: intermediate

.seealso: TaoView(), TaoLinearSolve()

.keywords: linear solver
@*/
int TaoViewLinearSolver(TAO_SOLVER solver){
  int info;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  if (solver->ksp){
    info = solver->ksp->View();CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetLinearSolverOptions"
/*@C
   TaoSetLinearSolverOptions - Set options for the linear solver

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Level: developer

.keywords: line search, options

.seealso: TaoGetLinearSolver()
@*/
int TaoSetLinearSolverOptions(TAO_SOLVER tao){
  int info;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (tao->ksp){
    info = tao->ksp->SetOptions();CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

