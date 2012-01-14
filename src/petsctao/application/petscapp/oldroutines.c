#include "tao.h"   /*I  "tao.h"  I*/


class TaoPetscApplication;
int TaoAppGetTaoPetscApp(TAO_APPLICATION, TaoPetscApplication**);

typedef struct {
  /* Function Gradient Evaluation over single element  */
  int  (*computef)(TAO_SOLVER, Vec, double*, void*); 
  int  (*computeg)(TAO_SOLVER, Vec, Vec, void*); 
  int  (*computefg)(TAO_SOLVER, Vec, double*, Vec, void*); 
  int  (*computeh)(TAO_SOLVER, Vec, Mat*, Mat*, MatStructure*, void*); 
  int  (*computer)(TAO_SOLVER, Vec, Vec, void*); 
  int  (*computej)(TAO_SOLVER, Vec, Mat*, void*); 
  TAO_SOLVER tao;
  void *fctx;
  void *fgctx;
  void *gctx;
  void *hctx;
  void *rctx;
  void *jctx;
} TaoPETScOLD;


#undef __FUNCT__  
#define __FUNCT__ "TaoPetscAppF"
static int TaoPetscAppF(TAO_APPLICATION taoapp , Vec X , double *f, void*ctx){
  int info;
  TaoPETScOLD* mctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  info=(*mctx->computef)(mctx->tao,X,f,mctx->fctx); CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscAppG"
static int TaoPetscAppG(TAO_APPLICATION taoapp , Vec X , Vec G, void*ctx){
  int info;
  TaoPETScOLD* mctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  info=(*mctx->computeg)(mctx->tao,X,G,mctx->gctx); CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscAppFG"
static int TaoPetscAppFG(TAO_APPLICATION taoapp , Vec X , double *f, Vec G, void* ctx){
  int info;
  TaoPETScOLD* mctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  info=(*mctx->computefg)(mctx->tao,X,f,G,mctx->fgctx); CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscAppH"
static int TaoPetscAppH(TAO_APPLICATION taoapp , Vec X , Mat *H, Mat *HP, MatStructure *flag, void* ctx){
  int info;
  TaoPETScOLD* mctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  info=(*mctx->computeh)(mctx->tao,X,H,HP,flag,mctx->hctx); CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscAppR"
static int TaoPetscAppR(TAO_APPLICATION taoapp , Vec X , Vec R, void*ctx){
  int info;
  TaoPETScOLD* mctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  info=(*mctx->computer)(mctx->tao,X,R,mctx->rctx); CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscAppJ"
static int TaoPetscAppJ(TAO_APPLICATION taoapp , Vec X , Mat *J, Mat *JP, MatStructure *flag, void* ctx){
  int info;
  TaoPETScOLD* mctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  info=(*mctx->computej)(mctx->tao,X,J,mctx->jctx); CHKERRQ(info);
  PetscFunctionReturn(0);  
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetPetscFunction"
/* @C
   TaoSetPetscFunction - Sets a routine that evaluates the function at
the specified point.

   Collective on TAO_SOLVER

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  x - variable vector that stores the solution
.  func - function evaluation routine
-  ctx - [optional] user-defined context for private data for the 
         function and gradient evaluation routine (may be TAO_NULL)

   Note: This function is no longer supported. Use TaoSetObjectiveFunction() and TaoSetVariableVec() instead.

   Level: intermediate

.seealso:  TaoAppSetObjectiveRoutine(), TaoAppSetInitialSolutionVec()
@ */
int TaoSetPetscFunction(TAO_APPLICATION taoapp, Vec X, int (*func)(TAO_SOLVER,Vec,double*,void*),void *ctx){
  int info;
  TaoPETScOLD* fgctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  PetscNew(TaoPETScOLD,&fgctx);
  fgctx->computef = func;
  fgctx->fctx=ctx;
  fgctx->tao = 0;
  info = TaoAppSetInitialSolutionVec(taoapp,X); CHKERRQ(info);
  info = TaoAppSetObjectiveRoutine(taoapp,TaoPetscAppF,(void*)fgctx); CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(taoapp,TaoApplicationFreeMemory, (void*)fgctx); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetPetscGradient"
/* @C
   TaoSetPetscGradient - Sets the gradient evaluation routine and gradient
   vector for use by the TAO_SOLVER routines.

   Collective on TAO_SOLVER

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  g - vector to store gradient
.  grad - gradient evaluation routine
-  ctx - [optional] user-defined function context 

   Note:
   This routine is no longer supported.

.seealso: TaoAppSetGradientRoutine()

@ */
int TaoSetPetscGradient(TAO_APPLICATION taoapp, Vec g, int (*grad)(TAO_SOLVER,Vec,Vec,void*),void *ctx){
  int info;
  TaoPETScOLD* fgctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  PetscNew(TaoPETScOLD,&fgctx);
  fgctx->computeg = grad;
  fgctx->gctx=ctx;
  fgctx->tao = 0;
  info = TaoAppSetGradientRoutine(taoapp,TaoPetscAppG,(void*)fgctx); CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(taoapp,TaoApplicationFreeMemory, (void*)fgctx); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetPetscFunctionGradient"
/* @C
   TaoSetPetscFunctionGradient - Sets a routine for function and gradient evaluation.

   Collective on TAO_SOLVER

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  x - vector to store solution
.  g - vector to store gradient
.  funcgrad - routine for evaluating the function and gradient
-  ctx - optional user-defined context for private data for the 
         function and gradient evaluation routine (may be TAO_NULL)

   Level: intermediate

   Note: 
   This routine is no longer supported.

.seealso: TaoAppSetObjectiveAndGradientRoutine(), TaoSetVariableVec()

@ */
int TaoSetPetscFunctionGradient(TAO_APPLICATION taoapp, Vec x, Vec g, int (*funcgrad)(TAO_SOLVER,Vec,double*,Vec, void*),void *ctx){
  int info;
  TaoPETScOLD* fgctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  PetscNew(TaoPETScOLD,&fgctx);
  fgctx->computefg = funcgrad;
  fgctx->fgctx=ctx;
  fgctx->tao = 0;
  info = TaoAppSetObjectiveAndGradientRoutine(taoapp,TaoPetscAppFG,(void*)fgctx); CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(taoapp,TaoApplicationFreeMemory, (void*)fgctx); CHKERRQ(info);
  info = TaoAppSetInitialSolutionVec(taoapp,x); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetPetscHessian"
/* @C
   TaoSetPetscHessian - Sets the function to compute the Hessian as well as the
   location to store the matrix.

   Collective on TAO_SOLVER and Mat

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  H - Hessian matrix
.  Hpre - preconditioner matrix (usually same as the Hessian)
.  hess - Hessian evaluation routine
-  ctx - [optional] user-defined context for private data for the 
         Hessian evaluation routine (may be TAO_NULL)

   Level: intermediate

   Note: This routine is no longer supported.

.seealso: TaoAppSetHessianRoutine(), TaoAppSetHessianMat()
@ */
int TaoSetPetscHessian(TAO_APPLICATION taoapp, Mat H, Mat Hpre, int (*hess)(TAO_SOLVER,Vec,Mat*,Mat*,MatStructure*,void*),void *ctx){
  int info;
  TaoPETScOLD* fgctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  PetscNew(TaoPETScOLD,&fgctx);
  fgctx->computeh = hess;
  fgctx->hctx=ctx;
  fgctx->tao = 0;
  info = TaoAppSetHessianRoutine(taoapp,TaoPetscAppH,(void*)fgctx);CHKERRQ(info);
  info = TaoAppSetHessianMat(taoapp,H,Hpre);CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(taoapp,TaoApplicationFreeMemory, (void*)fgctx); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetPetscConstraintsFunction"
/* @C
   TaoSetPetscConstraintsFunction - Sets the routine that evaluates
equality constraint functions.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  r - vector to constrainf function values
.  func - the constraint function evaluation routine
-  ctx - [optional] user-defined function context 

   Calling sequence of func:
$     func (TAO_SOLVER solver,Vec x,Vec r,void *ctx);

+  solver - the TAO_SOLVER solver context
.  x - input vector
.  r - constraint vector
-  ctx - user-defined function gradient context set from TaoSetPetscFunction()

   Level: intermediate

.keywords: TAO_SOLVER, set, function

.seealso: TaoGetGradient(), TaoSetPetscHessian()

@ */
int TaoSetPetscConstraintsFunction(TAO_APPLICATION taoapp, Vec r, int (*func)(TAO_SOLVER,Vec,Vec,void*),void 
				   *ctx){
  int info;
  TaoPETScOLD* fgctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  PetscNew(TaoPETScOLD,&fgctx);
  fgctx->computer = func;
  fgctx->rctx=ctx;
  fgctx->tao = 0;
  info = TaoAppSetConstraintRoutine(taoapp,TaoPetscAppR,(void*)fgctx); CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(taoapp,TaoApplicationFreeMemory, (void*)fgctx); CHKERRQ(info);
  info = TaoAppSetFunctionVec(taoapp,r); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetPetscJacobian"
/* @C
   TaoSetPetscJacobian - Sets the function to compute the Jacobian of
the equality constraint function as well as the
   location to store the matrix.

   Collective on TAO_APPLICATION and Mat

   Input Parameters:
+  solver - the TAO_APPLICATION context
.  J - Jacobian matrix
.  jac - Jacobian evaluation routine
-  ctx - [optional] user-defined context for private data for the 
         Hessian evaluation routine (may be TAO_NULL)

   Calling sequence of func:
$    jac (TAO_SOLVER solver,Vec x,Mat *J,void *ctx);

+  solver - the TAO_SOLVER solver context
.  x - input vector
.  J - Jacobian matrix
-  ctx - [optional] user-defined Hessian context

   Notes: 

   The function jac() takes Mat * as the matrix arguments rather than Mat.  
   This allows the Jacobian evaluation routine to replace J with a 
   completely new new matrix structure (not just different matrix elements)
   when appropriate, for instance, if the nonzero structure is changing
   throughout the global iterations.

   Level: intermediate

.keywords: TAO_SOLVER, Jacobian

.seealso: TaoSetPetscFunction(), TaoSetPetscConstraintsFunction()
@ */
int TaoSetPetscJacobian(TAO_APPLICATION taoapp,Mat J,int (*jac)(TAO_SOLVER,Vec,Mat*,void*),void *ctx){
  int info;
  TaoPETScOLD* fgctx= (TaoPETScOLD*)ctx;
  PetscFunctionBegin;
  PetscNew(TaoPETScOLD,&fgctx);
  fgctx->computej = jac;
  fgctx->jctx=ctx;
  fgctx->tao = 0;
  info = TaoAppSetJacobianRoutine(taoapp,TaoPetscAppJ,(void*)fgctx);CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(taoapp,TaoApplicationFreeMemory, (void*)fgctx); CHKERRQ(info);
  info = TaoAppSetJacobianMat(taoapp,J,J); CHKERRQ(info);
  PetscFunctionReturn(0);
}




#undef __FUNCT__  
#define __FUNCT__ "TaoSetFromOptions_Old"
int TaoSetFromOptions(TAO_APPLICATION taoapp){
  int info;
  PetscFunctionBegin;
  info=TaoAppSetFromOptions(taoapp); CHKERRQ(info);
  PetscFunctionReturn(0);
}

class TaoApplication;
extern int TaoSetApplication(TAO_SOLVER, TaoApplication*);

#undef __FUNCT__  
#define __FUNCT__ "TaoSetApplication_Old"
int TaoSetApplication(TAO_SOLVER tao, TAO_APPLICATION taoapp){
  int info;
  TaoPetscApplication* tpapp;
  PetscFunctionBegin;
  info=TaoAppGetTaoPetscApp(taoapp, &tpapp); CHKERRQ(info);
  info=TaoSetApplication(tao,(TaoApplication*)(tpapp)); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoApplicationDestroy_OLD"
int TaoApplicationDestroy(TAO_APPLICATION taoapp){
  int info;
  PetscFunctionBegin;
  info=TaoAppDestroy(taoapp); CHKERRQ(info);
  PetscFunctionReturn(0);
}
