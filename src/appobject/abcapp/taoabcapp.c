#include "src/tao_impl.h"     /*I  "tao_solver.h"  I*/
#include "taoabcapp.h"
#include "src/vector/tvecdouble.h"

#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::TaoABCApplication"
TaoABCApplication::TaoABCApplication(){
  this->taox=0;
  return;
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::~TaoABCApplication"
TaoABCApplication::~TaoABCApplication(){
  int info;
  info = TaoVecDestroy(this->taox);
  if (info) return;
  return;
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::SetDimension"
/* @C
   SetDimension - Set the number of variables in this application.

   Input Parameters:
.  n - the number of variables

   Note:  This method should be called only once and be called before 
   the application is set to the solver.

   Note: This routine creates a variables vector.  Applications derived
   from this object may create a TaoVec object and set the 'taox' field
   another way.

   Level: beginner

@ */
int TaoABCApplication::SetNumberOfVariables(TaoInt n){
  TaoFunctionBegin; 
  this->taox = new TaoVecDoubleArray(n);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::GetSolutionAndGradient"
/* @C
   GetSolutionAndGradient - Get the solution and and gradient arrays.

   Output Parameters:

+  n - the length of the arrays, which equals the number of variables
.  x - variables vector
-  g - gradient vector

   Level: beginner

@ */
int TaoABCApplication::GetSolution(double* &x, TaoInt &n){
  TaoInt nn;
  int info;
  double *xx;
  TaoFunctionBegin;
  info = this->taox->GetDoubles(&xx,&nn);CHKERRQ(info);
  x=xx; n=nn;
  info = this->taox->RestoreDoubles(&xx,&nn);  CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::GetVariableVector"
int TaoABCApplication::GetVariableVector(TaoVec **xx){

  TaoFunctionBegin;
  *xx= this->taox;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::EvaluateFunction"
int TaoABCApplication::EvaluateObjectiveFunction(TaoVec *xx, double *ff){
  TaoFunctionBegin;
  TaoFunctionReturn(1);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::EvaluateGradient"
int TaoABCApplication::EvaluateGradient(TaoVec *tx, TaoVec *tg){
  TaoInt   n;
  int info;
  double *xptr,*gptr;
  double ff;
  TaoVecDoubleArray *xx=(TaoVecDoubleArray*)tx, *gg=(TaoVecDoubleArray*)tg;
  TaoFunctionBegin;
  info = xx->GetDoubles(&xptr,&n);CHKERRQ(info);
  info = gg->GetDoubles(&gptr,&n);CHKERRQ(info);
  info = this->ComputeObjectiveAndGradient(xptr,n,&ff,gptr);CHKERRQ(info);
  info = gg->RestoreDoubles(&gptr,&n);CHKERRQ(info);
  info = xx->RestoreDoubles(&xptr,&n);  CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::EvaluateObjectiveAndGradient"
int TaoABCApplication::EvaluateObjectiveAndGradient(TaoVec *tx, double *ff, TaoVec *tg){
  TaoInt     n;
  int info;
  double *xptr, *gptr;
  TaoVecDoubleArray *xx=(TaoVecDoubleArray*)tx, *gg=(TaoVecDoubleArray*)tg;

  TaoFunctionBegin;
  info = xx->GetDoubles(&xptr,&n);CHKERRQ(info);
  info = gg->GetDoubles(&gptr,&n);CHKERRQ(info);
  info = this->ComputeObjectiveAndGradient(xptr,n,ff,gptr);CHKERRQ(info);
  info = gg->RestoreDoubles(&gptr,&n);CHKERRQ(info);
  info = xx->RestoreDoubles(&xptr,&n);  CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::ComputeObjectiveAndGradient"
/*@C
   ComputeObjectiveAndGradient - Compute the objective function and its gradient.

   Input Parameters:
.  x - an array with the point at which the function and gradient should be evalutated.
.  g - an array that will contain the gradient vector.
+  n - the length of the arrays, which equals the number of variables

   Output Parameters:

.  f - function value
.  g - gradient vector

   Level: beginner

@*/
// virtual int TaoABCApplication::ComputeObjectiveAndGradient(double*x,TaoInt n,double*f,double *g){}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::StartingPoint"
/*@C
   StartingPoint - Define the starting point for the solver.

   Collective on TAO_SOLVER

   Input Parameters:
+  x - vector to store initial variables
-  n - length of variable vector array

   Note: This virtual method should be implemented in the application
   if a starting point other than 0 is desired.

   Level: beginner

@*/
int TaoABCApplication::StartingPoint(double *x, TaoInt n){
  TaoInt  i;
  TaoFunctionBegin;
  for (i=0;i<n;i++) x[i]=0.0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::InitializeVariables"
int TaoABCApplication::InitializeVariables(TaoVec *tx){
  TaoInt n;
  int info;
  double *xptr;
  TaoVecDoubleArray *xx=(TaoVecDoubleArray*)tx; 
  TaoFunctionBegin; 
  info = xx->GetDoubles(&xptr,&n);CHKERRQ(info);
  info = this->StartingPoint(xptr,n);CHKERRQ(info);
  info = xx->RestoreDoubles(&xptr,&n);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::EvaluateVariableBounds"
int TaoABCApplication::EvaluateVariableBounds(TaoVec *txl, TaoVec *txu){
  TaoInt n;
  int info;
  double *xl,*xu;
  TaoVecDoubleArray *xxll=(TaoVecDoubleArray*)txl, *xxuu=(TaoVecDoubleArray*)txu;
  TaoFunctionBegin;
  
  info = xxll->GetDoubles(&xl,&n);CHKERRQ(info);
  info = xxuu->GetDoubles(&xu,&n);CHKERRQ(info);
  info = this->ComputeVariableBounds(xl, xu,n);CHKERRQ(info);
  info = xxll->RestoreDoubles(&xu,&n);CHKERRQ(info);
  info = xxuu->RestoreDoubles(&xl,&n);  CHKERRQ(info);

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCApplication::ComputeVariableBounds"
/*@C
   StartingPoint - Define the lower and upper bounds on the variables.

   Collective on TAO_SOLVER

   Input Parameters:
+  xl - array to store lower bounds
.  xu - array to store upper bounds
-  n - length of these arrays, which equals the number of variables

   Note:  This virtual method should be implemented by the application
   if there are lower or upper bounds on the variables.

   Level: beginner

@*/
int TaoABCApplication::ComputeVariableBounds(double *xl, double *xu, TaoInt n){
  TaoInt  i;
  TaoFunctionBegin;
  for (i=0;i<n;i++){ 
    xl[i]=TAO_NINFINITY;
    xu[i]=TAO_INFINITY;
  }
  TaoFunctionReturn(0);
}
