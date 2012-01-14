#include "src/tao_impl.h"     /*I  "tao_solver.h"  I*/
#include "taofloatapp.h"
#include "src/vector/tvecsingle.h"

#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::TaoABCFloatApplication"
TaoABCFloatApplication::TaoABCFloatApplication(){
  this->taox=0;
  return;
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::~TaoABCFloatApplication"
TaoABCFloatApplication::~TaoABCFloatApplication(){
  int info;
  info = TaoVecDestroy(this->taox);
  if (info) return;
  return;
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::SetDimension"
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

@*/
int TaoABCFloatApplication::SetNumberOfVariables(TaoInt n){
  TaoFunctionBegin; 
  this->taox = new TaoVecFloatArray(n);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::GetSolutionAndGradient"
/* @C
   GetSolutionAndGradient - Get the solution and and gradient arrays.

   Output Parameters:

+  n - the length of the arrays, which equals the number of variables
.  x - variables vector
-  g - gradient vector

   Level: beginner

@ */
int TaoABCFloatApplication::GetSolutionAndGradient(float* &x, float* &g, TaoInt &n){
  TaoInt nn;
  int info;
  float *xx;
  TaoFunctionBegin;
  info = this->taox->GetFloats(&xx,&nn);CHKERRQ(info);
  x=xx; n=nn;
  g=0;
  info = this->taox->RestoreFloats(&xx,&nn);  CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::GetVariableVector"
int TaoABCFloatApplication::GetVariableVector(TaoVec **xx){

  TaoFunctionBegin;
  *xx= this->taox;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::EvaluateFunction"
int TaoABCFloatApplication::EvaluateObjectiveFunction(TaoVec *xx, double *ff){
  TaoFunctionBegin;
  TaoFunctionReturn(1);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::EvaluateGradient"
int TaoABCFloatApplication::EvaluateGradient(TaoVec *tx, TaoVec *tg){
  TaoInt   n;
  int info;
  float *xptr,*gptr;
  float ff;
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tx);
  TaoVecFloatArray* gg =  (TaoVecFloatArray*)(tg);
  TaoFunctionBegin;
  info = xx->GetFloats(&xptr,&n);CHKERRQ(info);
  info = gg->GetFloats(&gptr,&n);CHKERRQ(info);
  info = this->ComputeObjectiveAndGradient(xptr,n,&ff,gptr);CHKERRQ(info);
  info = gg->RestoreFloats(&gptr,&n);CHKERRQ(info);
  info = xx->RestoreFloats(&xptr,&n);  CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::EvaluateObjectiveAndGradient"
int TaoABCFloatApplication::EvaluateObjectiveAndGradient(TaoVec *tx, double *ff, TaoVec *tg){
  TaoInt     n;
  int info;
  float *xptr, *gptr;
  float f;
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tx);
  TaoVecFloatArray* gg =  (TaoVecFloatArray*)(tg);

  TaoFunctionBegin;
  info = xx->GetFloats(&xptr,&n);CHKERRQ(info);
  info = gg->GetFloats(&gptr,&n);CHKERRQ(info);
  info = this->ComputeObjectiveAndGradient(xptr,n,&f,gptr);CHKERRQ(info);
  info = gg->RestoreFloats(&gptr,&n);CHKERRQ(info);
  info = xx->RestoreFloats(&xptr,&n);  CHKERRQ(info);
  *ff=(double)f;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::ComputeObjectiveAndGradient"
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
// virtual int TaoABCFloatApplication::ComputeObjectiveAndGradient(float*x,TaoInt n,float*f,float *g){}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::StartingPoint"
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
int TaoABCFloatApplication::StartingPoint(float *x, TaoInt n){
  TaoInt  i;
  TaoFunctionBegin;
  for (i=0;i<n;i++) x[i]=0.0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::InitializeVariables"
int TaoABCFloatApplication::InitializeVariables(TaoVec *tx){
  TaoInt n;
  int info;
  float *xptr;
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tx);
  TaoFunctionBegin; 
  info = xx->GetFloats(&xptr,&n);CHKERRQ(info);
  info = this->StartingPoint(xptr,n);CHKERRQ(info);
  info = xx->RestoreFloats(&xptr,&n);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::EvlauateVariableBounds"
int TaoABCFloatApplication::EvaluateVariableBounds(TaoVec *txl, TaoVec *txu){
  TaoInt n;
  int info;
  float *xl,*xu;
  TaoVecFloatArray* xxll =  (TaoVecFloatArray*)(txl);
  TaoVecFloatArray* xxuu =  (TaoVecFloatArray*)(txu);

  TaoFunctionBegin;
    
  info = xxll->GetFloats(&xl,&n);CHKERRQ(info);
  info = xxuu->GetFloats(&xu,&n);CHKERRQ(info);
  info = this->ComputeVariableBounds(xl, xu,n);CHKERRQ(info);
  info = xxuu->RestoreFloats(&xu,&n);CHKERRQ(info);
  info = xxll->RestoreFloats(&xl,&n);  CHKERRQ(info);
  
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoABCFloatApplication::ComputeVariableBounds"
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
int TaoABCFloatApplication::ComputeVariableBounds(float *xl, float *xu, TaoInt n){
  TaoInt  i;
  TaoFunctionBegin;
  for (i=0;i<n;i++){ 
    xl[i]=TAO_NINFINITY;
    xu[i]=TAO_INFINITY;
  }
  TaoFunctionReturn(0);
}
