#include "taoappobject.h"     /*I "tao_solver.h"  I*/
#include "tao_general.h"
#include "taovec.h"


#undef __FUNCT__
#define __FUNCT__ "TaoApplication::EvaluateObjectiveFunction"
/*@C
   EvaluateObjectiveFunction - Evaluate the objective function
   at the point x.

   Collective on TAO_SOLVER

   Input Parameters:
.  xx - the point at which to evaluate the obective function

   Output Parameters:
.  ff - objective function value

.seealso TaoSetPetscGradient(), EvaluateObjectiveFunction(), EvaluateObjectiveAndGradientFunction()

   Level: intermediate

.keywords: application, context

@*/
int TaoApplication::EvaluateObjectiveFunction(TaoVec *xx, double *ff){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  //  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::EvaluateGradient"
/*@C
   EvaluateGradient - Evaluate the gradient of the
   objective function at the point x.

   Collective on TAO_SOLVER

   Input Parameters:
.  xx - the point at which to evaluate the gradient

   Output Parameters:
.  gg - the gradient vector

.seealso TaoSetPetscGradient(), EvaluateObjectiveFunction(), EvaluateObjectiveAndGradientFunction()

   Level: intermediate

.keywords: application, context

@*/
int TaoApplication::EvaluateGradient(TaoVec *xx, TaoVec *gg){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  //  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoApplication::EvaluateObjectiveAndGradient"
/*@C
   EvaluateObjectiveAndGradient - Evaluate the objective function
   and its gradient at the point x.

   Collective on TAO_SOLVER

   Input Parameters:
.  xx - the point at which to evaluate the gradient

   Output Parameters:
+  ff - the objective value
-  gg - the gradient vector

.seealso TaoSetPetscFunctionGradient(), EvaluateGradient(), EvaluateObjectiveFunction()

   Level: intermediate

.keywords: application, gradient

@*/
int TaoApplication::EvaluateObjectiveAndGradient(TaoVec *xx, double *ff, TaoVec *gg){
  int info;
  TaoFunctionBegin;
  info=this->EvaluateObjectiveFunction(xx,ff);CHKERRQ(info);
  info=this->EvaluateGradient(xx,gg);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::EvaluateHessian"
/*@C
   EvaluateHessian - Evaluate the Hessian of the objective function
   at the point x.

   Collective on TAO_SOLVER

   Input Parameters:
.  xx - the point at which to evaluate the gradient

   Output Parameters:
.  HH - the Hessian matrix

.seealso TaoSetPetscHessian()
   Level: intermediate

.keywords: application, hessian

@*/
int TaoApplication::EvaluateHessian(TaoVec *xx, TaoMat *HH){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  //  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::HessianSolve"
/*@C
   HessianSolve - Solve a linear system involving the Hessian
matrix, or apply an inverse Hessian operator to a vector.

   Collective on TAO_SOLVER

   Input Parameters:
.  rhs - the right-hand side

   Output Parameters:
+  vv - the solution
-  success - flag states whether solution was found.

.seealso TaoSetMethod(), TaoLMVMSetSize(), EvaluateHessian()

   Level: intermediate

.keywords: application, Hessian, LMVM

@*/
int TaoApplication::HessianSolve(TaoVec *rhs, TaoVec*vv, TaoTruth *success){
  int info;
  TaoFunctionBegin;
  info=vv->CopyFrom(rhs);CHKERRQ(info);
  *success=TAO_TRUE;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoApplication::EvaluateConstraints"
/*@C
   EvaluateConstraints - Evaluate the constraint functions
   at the point x.

   Collective on TAO_SOLVER

   Input Parameters:
.  xx - the point at which to evaluate the gradient

   Output Parameters:
.  RR - the constraint vector

.seealso TaoSetPetscConstraintsFunction()
   Level: intermediate

.keywords: application, constraints

@*/
int TaoApplication::EvaluateConstraints(TaoVec *xx, TaoVec *RR){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  //  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::EvaluateJacobian"
/*@C
   EvaluateJacobian - Evaluate the Jacobian of the constraint functions
   at the point x.

   Collective on TAO_SOLVER

   Input Parameters:
.  xx - the point at which to evaluate the gradient

   Output Parameters:
.  JJ - the Jacobian matrix

.seealso TaoSetPetscJacobian()

   Level: intermediate

.keywords: application, constraints

@*/
int TaoApplication::EvaluateJacobian(TaoVec *xx, TaoMat *JJ){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  //  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::InitializeVariables"
/*@C
   InitializeVariables - Initialize the variables for the optimization
   solver.  Set an initial point.


   Collective on TAO_SOLVER

   Input Parameters:
.  xx - the variable vector

   Level: intermediate

.keywords: application, constraints

@*/
int TaoApplication::InitializeVariables(TaoVec *xx){
  int info;
  TaoFunctionBegin;
  info=xx->SetToZero();CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::GetVariableVector"
/*@C
   GetVariableVector - Sets the pointer to a vector that will
   be used to store the variables.

   Collective on TAO_SOLVER

   Output Parameters:
.  xx - vector to the variable vector

   Level: intermediate

.keywords: application, variables

@*/
int TaoApplication::GetVariableVector(TaoVec **xx){
  TaoFunctionBegin;
  *xx=0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoApplication::EvaluateVariableBounds"
/*@C
   EvaluateVariableBounds - Set these vectors equal to the lower
   and upper bounds on the variables.

   Collective on TAO_SOLVER

   Input Parameters:
+  xxll - vector vector of lower bounds
-  xxuu - vector vector of upper bounds

   Level: intermediate

.keywords: application, bounds

@*/
int TaoApplication::EvaluateVariableBounds(TaoVec *xxll, TaoVec *xxuu)
{
  int info;

  TaoFunctionBegin;

  if (xxll) {
    info = xxll->SetToConstant(-TAO_INFINITY); CHKERRQ(info);
  }

  if (xxuu) {
    info = xxuu->SetToConstant(TAO_INFINITY); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoApplication::GetHessianMatrix"
/*@C
   GetHessianMatrix - Sets the pointer to a matrix that will
   be used to store the Hessian of the objective function.

   Collective on TAO_SOLVER

   Output Parameters:
.  HH - vector to the Hessian

   Level: intermediate

.keywords: application, gradient

@*/
int TaoApplication::GetHessianMatrix(TaoMat **HH){
  TaoFunctionBegin;
  *HH=0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::GetJacobianMatrix"
/*@C
   GetJacobianMatrix - Sets the pointer to a matrix that will
   be used to store the Jacobian

   Collective on TAO_SOLVER

   Output Parameters:
.  JJ - matrix to the Jacobian

   Level: intermediate

.keywords: application, gradient

@*/
int TaoApplication::GetJacobianMatrix(TaoMat **JJ){
  TaoFunctionBegin;
  *JJ=0;
  TaoFunctionReturn(0);
}

int TaoApplication::GetConstraintVector(TaoVec **rr){
  TaoFunctionBegin;
  *rr=0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoApplication::GetInequalityConstraints"
int TaoApplication::GetInequalityConstraints(TaoVec**ll, TaoMat **AA, TaoVec **uu){
  TaoFunctionBegin;
  *ll=0;
  *AA=0;
  *uu=0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::Monitor"
/*@C
   Monitor - Monitor the current solution

   Collective on TAO_SOLVER

   Level: intermediate

   Note: This routine is called after each iteration

.keywords: application, monitor

@*/
int TaoApplication::Monitor(){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApplication::Monitor2"
/*@C
   Monitor2 - 

   Collective on TAO_SOLVER

   Input Parameters:
.  xx - the current solution
.  gl - the gradient of the Lagrangian function
.  dx - the step direction

   Output Parameters:
.  *term - TAO_TRUE if the solver should stop iterating and TAO_FALSE if the solver should continue;

.seealso TaoSetPetscGradient(),

   Level: intermediate

.keywords: application, context

@*/
int TaoApplication::Monitor2(TaoVec *xx, TaoVec *gg, TaoVec *ff, TaoTruth *term){
  TaoFunctionBegin;
  *term=TAO_FALSE;
  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "TaoApplication::GetLinearSolver"
/*@C
   GetLinearSolver - Create and a linear solver used to solve this matrix

   Collective on TAO_SOLVER

   input Parameters:
+  H - the operator to be used to solve
-  flag - indicates properties solver must have

   Output Parameters:
.  tksp - the linear solver

   Types of linear solvers are:
+  100 - Solver for any square matrix
.  110 - GMRES
.  200 - Any symmetric matrix
.  210 - MINRES
.  220 - CG with Trust  Region
.  300 - Any symmetric positive definite
-  310 - Conjugate Gradient

   Level: intermediate

.keywords: application, linear solver

@*/
int TaoApplication::GetLinearSolver(TaoMat *H, int flag, TaoLinearSolver **tksp){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
}

#undef __FUNCT__  
#define __FUNCT__ "TaoDestroyApplication"
/*@
   TaoDestroyApplication - Destroys the TAO application

   Collective on TAO_APPLICATION

   Input Parameters:
.  myapp - pointer to the TaoApplication object (TaoApplication*),


   Level: developer

.keywords: application

@*/
int TaoDestroyApplication(TaoApplication *myapp){
  TaoFunctionBegin;
  delete myapp;
  TaoFunctionReturn(0);
}


