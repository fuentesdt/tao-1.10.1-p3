#include "tao_general.h"
#include "tao_app_impl.h"     /*I  "tao.h"  I*/
#include "src/petsctao/include/taopetsc.h"

extern int Tao_ObjectiveEval, Tao_GradientEval, Tao_HessianEval;
extern int TAO_APP_COOKIE;

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetObjectiveRoutine"
/*@C
   TaoAppSetObjectiveRoutine - Sets a routine that evaluates the function at
the specified point.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  func - function evaluation routine
-  ctx - [optional] user-defined context for private data for the 
         function and gradient evaluation routine (may be TAO_NULL)

   Calling sequence of func:
$     func (TAO_APPLICATION taoapp, Vec x,double *f,void *ctx);

+  taoapp - the TAO_APPLICATION context
.  x - input vector
.  f - objective function value
-  ctx - [optional] user-defined function context 

   Note:  
   Most applications do not need this routine.  The routine
   TaoAppSetObjectiveFunctionGradient() is sufficient.

   Level: intermediate

.keywords: TAO_APPLICATION, set, minimization, function

.seealso:  TaoAppSetHessianRoutine(), TaoAppSetObjectiveAndGradientRoutine(), TaoAppSetInitialSolutionVec()
@*/
int TaoAppSetObjectiveRoutine(TAO_APPLICATION taoapp, int (*func)(TAO_APPLICATION,Vec,double*,void*),void *ctx){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->computeumfunction=func;
  taoapp->usrfctx=ctx;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppComputeObjective"
/*@
   TaoAppComputeObjective - Compute the objective function that has been
   set with TaoAppSetObjectiveRoutine().

   Collective on TAO_APPLICATION

   Input Parameters:
+  taopp - the TAO_APPLICATION context
-  X - the point where the objective should be evaluated

   Output Parameter:
.  f - function value

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeGradient(), TaoAppComputeObjectiveAndGradient()
@*/
int TaoAppComputeObjective(TAO_APPLICATION taoapp, Vec X, double *f){
  int     info;
  Vec G;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  PetscValidHeaderSpecific(X,VEC_COOKIE,2);
  PetscStackPush("TAO User minimization function");
  info = PetscLogEventBegin(Tao_ObjectiveEval,taoapp,X,0,0);
  if (taoapp->computeumfunction){
    info = (*taoapp->computeumfunction)(taoapp,X,f,taoapp->usrfctx); CHKERRQ(info);
  } else if (taoapp->computefunctiongradient){
    info = VecDuplicate(X,&G);CHKERRQ(info);
    info = (*taoapp->computefunctiongradient)(taoapp,X,f,G,taoapp->usrfgctx);
    CHKERRQ(info);
    info=VecDestroy(G);CHKERRQ(info);
    taoapp->ngeval++;
  } else {
    SETERRQ(1,"TAO ERROR: Must set Objective function");
  }
  taoapp->nfeval++;
  info = PetscLogEventEnd(Tao_ObjectiveEval,taoapp,X,0,0);
  PetscStackPop;
  info = PetscInfo1(taoapp,"TAO Function evaluation: %14.12e\n",*f);CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetGradientRoutine"
/*@C
   TaoAppSetGradientRoutine - Sets the gradient evaluation routine
   for use by the TAO_APPLICATION routines.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  grad - gradient evaluation routine
-  ctx - [optional] user-defined function context 

   Calling sequence of func:
$     grad (TAO_APPLICATION taoapp,Vec x,Vec g,void *ctx);

+  taoapp - the TAO_APPLICATION  context
.  x - input vector
.  g - gradient vector
-  ctx - user-defined function gradient context set from TaoAppSetGradientRoutine()

   Level: intermediate

   Options Database Keys:
.  -tao_view_gradient - view the gradient after each evaluation using PETSC_VIEWER_STDOUT_WORLD

   Note:
   In most cases, the routine TaoAppSetObjectiveAndGradientRoutine() is more appropriate.  However,
   when using finite differences to compute the Hessian, setting this routine can  be
   beneficial.

.keywords: TAO_APPLICATION, set, gradient

.seealso: TaoAppGetGradientVec(), TaoAppSetObjectiveAndGradientRoutine(), TaoAppSetHessianRoutine()

@*/
int TaoAppSetGradientRoutine(TAO_APPLICATION taoapp, int (*grad)(TAO_APPLICATION,Vec,Vec,void*),void *ctx){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->computegradient=grad;
  taoapp->usrgctx=ctx;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppComputeGradient"
/*@
   TaoAppComputeGradient - Compute the gradient of the objective function using the
   routine set by TaoApplicationSetGradientRoutine().

   Collective on TAO_APPLICATION

   Input Parameters:
+  taopp - the TAO_APPLICATION context
-  X - the point where the objective should be evaluated

   Output Parameter:
.  f - function value

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeGradient(), TaoAppComputeObjectiveAndGradient()
@*/
int TaoAppComputeGradient(TAO_APPLICATION taoapp, Vec X, Vec G){
  int     info;
  double ff;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  PetscValidHeaderSpecific(G,VEC_COOKIE,3);
  PetscValidHeaderSpecific(X,VEC_COOKIE,2);
  PetscStackPush("TAO User Gradient Evaluation");
  info = PetscLogEventBegin(Tao_GradientEval,taoapp,X,G,0);
  if (taoapp->computegradient){
    info = (*taoapp->computegradient)(taoapp,X,G,taoapp->usrgctx);
    CHKERRQ(info);
  } else if ( taoapp->computefunctiongradient ) {
    info = (*taoapp->computefunctiongradient)(taoapp,X,&ff,G,taoapp->usrfgctx);
    CHKERRQ(info);
    taoapp->nfeval++;
  } else {
    SETERRQ(1,"TAO ERROR: Must set gradient");
  }
  taoapp->ngeval++;
  info = PetscLogEventEnd(Tao_GradientEval,taoapp,X,G,0);
  PetscStackPop;
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoAppComputeObjectiveAndGradient"
/*@
   TaoAppComputeObjectiveAndGradient - Compute the gradient of the objective function using the
   routine set by TaoApplicationSetGradientRoutine().

   Collective on TAO_APPLICATION

   Input Parameters:
+  taopp - the TAO_APPLICATION context
-  X - the point where the objective should be evaluated

   Output Parameter:
+  f - function value
-  G - the gradient vector.

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeGradient(), TaoAppSetObjectiveAndGradientRoutine()
@*/
int TaoAppComputeObjectiveAndGradient(TAO_APPLICATION taoapp, Vec X, double *f, Vec G){
  int     info;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,VEC_COOKIE,2);
  PetscValidHeaderSpecific(G,VEC_COOKIE,4);
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  PetscStackPush("TAO User Objective and gradient function");
  if (taoapp->computefunctiongradient){
    info = PetscLogEventBegin(Tao_GradientEval,taoapp,X,G,0);CHKERRQ(info);
    info = PetscLogEventBegin(Tao_ObjectiveEval,taoapp,X,G,0);
    info = (*taoapp->computefunctiongradient)(taoapp,X,f,G,taoapp->usrfgctx);
    CHKERRQ(info);
    info = PetscLogEventEnd(Tao_ObjectiveEval,taoapp,X,G,0);CHKERRQ(info);

    if (taoapp->computegradient) {
      info = PetscInfo(taoapp,"Computing gradient routine separately, ignoring gradient calculated in function/gradient evaluation routine."); CHKERRQ(info);
      info = (*taoapp->computegradient)(taoapp,X,G,taoapp->usrfgctx); 
      CHKERRQ(info);
    }
    info = PetscLogEventEnd(Tao_GradientEval,taoapp,X,G,0);CHKERRQ(info);
  } else if ( taoapp->computeumfunction && taoapp->computegradient ) {
    info = PetscLogEventBegin(Tao_ObjectiveEval,taoapp,X,G,0);
    info = (*taoapp->computeumfunction)(taoapp,X,f,taoapp->usrfctx);
    CHKERRQ(info);
    info = PetscLogEventEnd(Tao_ObjectiveEval,taoapp,X,G,0);CHKERRQ(info);
    info = PetscLogEventBegin(Tao_GradientEval,taoapp,X,G,0);CHKERRQ(info);
    info = (*taoapp->computegradient)(taoapp,X,G,taoapp->usrgctx);
    CHKERRQ(info);
    info = PetscLogEventEnd(Tao_GradientEval,taoapp,X,G,0);CHKERRQ(info);
  } else {
    SETERRQ(1,"TAO ERROR: Must set objective function and gradient.");
  }
  info = PetscInfo1(taoapp,"TAO Function evaluation: %14.12e\n",*f);CHKERRQ(info);
  taoapp->nfeval++;
  taoapp->ngeval++;

  PetscStackPop;
  PetscFunctionReturn(0);
}

 

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetObjectiveAndGradientRoutine"
/*@C
   TaoAppSetObjectiveAndGradientRoutine - Sets a routine for function and gradient evaluation.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  funcgrad - routine for evaluating the function and gradient
-  ctx - optional user-defined context for private data for the 
         function and gradient evaluation routine (may be TAO_NULL)

   Calling sequence of funcgrad:
$     funcgrad (TAO_APPLICATION tao,Vec x,double *f,Vec g,void *ctx);

+  tao - TAO_APPLICATION  context
.  x - input vector
.  f - function value
.  g - gradient vector
-  ctx - optional user-defined context 

   Notes:
   The user may call TaoAppSetObjectiveAndGradientRoutine() to set a routine
   that evaluates both the function and gradient.  Alternatively, the
   user may call both TaoAppSetObjectiveRoutine() and TaoAppSetGradientRoutine() to set
   separate routines for function and gradient evaluation.  

   Using a single routine to compute the function and gradient, as
   specified via TaoAppSetObjectiveAndGradientRoutine(), may enable better performance
   for applications in which many of the function and gradient computations
   are identical.

   Fortran Note:
   If your Fortran compiler does not recognize symbols over 31 characters in length, then
   use the identical routine with the shortened name TaoAppSetObjectiveAndGradientRo()


   Level: beginner

   Options Database Keys:
.   -tao_view_gradient - view the gradient after each iteration using PETSC_VIEWER_STDOUT_WORLD

.keywords: TAO_APPLICATION, set, objective, gradient

.seealso: TaoAppComputeObjectiveAndGradient()

@*/
int TaoAppSetObjectiveAndGradientRoutine(TAO_APPLICATION taoapp, int (*funcgrad)(TAO_APPLICATION,Vec,double*,Vec, void*),void *ctx){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->computefunctiongradient=funcgrad;
  taoapp->usrfgctx=ctx;
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoAppComputeHessian"
/*@C
   TaoAppComputeHessian - Compute the Hessian of the objective function using the
   routine set by TaoApplicationSetGradientRoutine().

   Collective on TAO_APPLICATION

   Input Parameters:
+  taopp - the TAO_APPLICATION context
.  X - the variable vector
.  H - the Hessian matrix
.  HP - the preconditioner for the Hessian matrix.
-  flag - flag used in KSPSetOperators()

   Output Parameter:
+  H - the Hessian matrix
.  HP - the preconditioner for the Hessian matrix.
-  flag - flag used in KSPSetOperators()

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeObjectiveAndGradient()
@*/
int TaoAppComputeHessian(TAO_APPLICATION taoapp, Vec X, Mat *HH, Mat *HHPre, MatStructure *flag){
  int     info;
  Mat H=*HH,HPre=*HHPre;
  MatStructure pflag=*flag;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,VEC_COOKIE,2);
  PetscValidHeaderSpecific(H,MAT_COOKIE,3);
  PetscValidHeaderSpecific(HPre,MAT_COOKIE,4);
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  info = PetscLogEventBegin(Tao_HessianEval,X,H,0,0);CHKERRQ(info);
  PetscStackPush("TAO User Hessian Routine");
  if (taoapp->computehessian){
    info = (*taoapp->computehessian)(taoapp,X,&H,&HPre,&pflag,taoapp->usrhctx);
    CHKERRQ(info);
  } else {
    SETERRQ(1,"TAO Error: No Hessian Routine Available");
  }
  *HH=H;*HHPre=HPre;*flag=pflag;
  taoapp->nheval++;
  info = PetscLogEventEnd(Tao_HessianEval,X,H,0,0);CHKERRQ(info);
  PetscStackPop;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetHessianRoutine"
/*@C
   TaoAppSetHessianRoutine - Sets the function to compute the Hessian as well as the
   location to store the matrix.

   Collective on TAO_APPLICATION and Mat

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  hess - Hessian evaluation routine
-  ctx - [optional] user-defined context for private data for the 
         Hessian evaluation routine (may be TAO_NULL)

   Calling sequence of hess:
$    hess (TAO_APPLICATION taoapp,Vec x,Mat *H,Mat *Hpre,int *flag,void *ctx);

+  taoapp - the TAO_APPLICATION  context
.  x - input vector
.  H - Hessian matrix
.  Hpre - preconditioner matrix, usually the same as A
.  flag - flag indicating information about the preconditioner matrix
   structure (see below)
-  ctx - [optional] user-defined Hessian context

   Options Database Keys:
.  -tao_view_hessian - view the hessian after each evaluation using PETSC_VIEWER_STDOUT_WORLD

   Notes: 

   The function hess() takes Mat * as the matrix arguments rather than Mat.  
   This allows the Hessian evaluation routine to replace A and/or B with a 
   completely new new matrix structure (not just different matrix elements)
   when appropriate, for instance, if the nonzero structure is changing
   throughout the global iterations.

   The flag can be used to eliminate unnecessary work in the preconditioner 
   during the repeated solution of linear systems of the same size.  The
   available options are
$    SAME_PRECONDITIONER -
$      B is identical during successive linear solves.
$      This option is intended for folks who are using
$      different Amat and Pmat matrices and want to reuse the
$      same preconditioner matrix.  For example, this option
$      saves work by not recomputing incomplete factorization
$      for ILU/ICC preconditioners.
$    SAME_NONZERO_PATTERN -
$      B has the same nonzero structure during
$      successive linear solves. 
$    DIFFERENT_NONZERO_PATTERN -
$      B does not have the same nonzero structure.

   Caution:
   If you specify SAME_NONZERO_PATTERN, the software believes your assertion
   and does not check the structure of the matrix.  If you erroneously
   claim that the structure is the same when it actually is not, the new
   preconditioner will not function correctly.  Thus, use this optimization
   feature carefully!

   If in doubt about whether your preconditioner matrix has changed
   structure or not, use the flag DIFFERENT_NONZERO_PATTERN.

   Level: beginner

.keywords: TAO_APPLICATION, Hessian

.seealso: TaoAppSetObjectiveAndGradientRoutine(), TaoAppSetHessianMat(), KSPSetOperators()
@*/
int TaoAppSetHessianRoutine(TAO_APPLICATION taoapp, int (*hess)(TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure*,void*),void *ctx){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->computehessian=hess;
  taoapp->usrhctx=ctx;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetHessianMat"
/*@
   TaoAppSetHessianMat - Sets the matrix representing the Hessian
   and the matrix used to precondition it.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  H - the matrix used for the Hessian
-  HP - the matrix used to precondition the Hessian matrix.

   Note:
   Usually H and HP are the same matrix

   Level: beginner

.seealso:  TaoAppSetHessianRoutine()
@*/
int TaoAppSetHessianMat(TAO_APPLICATION taoapp, Mat H, Mat HP ){
  int info;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);

  if (H) {
    PetscValidHeaderSpecific(H,MAT_COOKIE,2);
    PetscValidHeaderSpecific(HP,MAT_COOKIE,3);
    PetscObjectReference((PetscObject)H);
    PetscObjectReference((PetscObject)HP);
  }
  if (taoapp->H) {
    info=MatDestroy(taoapp->H);CHKERRQ(info); 
    info=MatDestroy(taoapp->HP);CHKERRQ(info); 
  }  

  taoapp->H=H;
  taoapp->HP=HP;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetHessianMat"
/*@
   TaoAppGetHessianMat - Sets the matrix representing the Hessian
   and the matrix used to precondition it.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  H - the matrix used for the Hessian
-  HP - the matrix used to precondition the Hessian matrix.

   Note:
   Usually H and HP are the same matrix

   Level: developer

.seealso:  TaoAppSetHessianMat()
@*/
int TaoAppGetHessianMat(TAO_APPLICATION taoapp, Mat *H, Mat *HP ){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (H) { *H=taoapp->H;}
  if (HP){ *HP=taoapp->HP;}
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetHessianSolveRoutine"
/*@C
   TaoAppSetHessianSolveRoutine - Sets the routine that solves a linear
system involving the Hessian operator, (or approximate Hessian).

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  ah - gradient evaluation routine
-  ctx - [optional] user-defined function context 

   Calling sequence of func:
$     ah (TAO_APPLICATION taoapp,Vec vin,Vec vout, PetscTruth* success, void *ctx);

+  taoapp - the TAO_APPLICATION  context
.  v - input vector
-  ctx - user-defined function gradient context set from TaoAppSetGradientRoutine()

   Level: intermediate

   Note:
   This routine is to be used only with the LMVM and BLMVM solvers.  These solvers do not use Hessian information, but can incorporate very approximate Hessian information into the routine.

.keywords: TAO_APPLICATION, set, hessian

.seealso: TaoAppSetHessianRoutine()

@*/
int TaoAppSetHessianSolveRoutine(TAO_APPLICATION taoapp, int (*ah)(TAO_APPLICATION,Vec,Vec,PetscTruth*,void*),void *ctx){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->hessiansolve=ah;
  taoapp->usrhhhctx=ctx;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppHessianSolve"
/*@
   TaoAppHessianSolve - Apply an
   inverse Hessian operator to the vector, or solve a linear
   system involving the Hessian.   It uses the
   routine set by TaoApplicationSetHessianSolveRoutine().

   Collective on TAO_APPLICATION

   Input Parameters:
+  taopp - the TAO_APPLICATION context
-  Vin - the vector to be applied to the approximate inverse Hessian operator.

   Output Parameter:
+  Vout - the inverse Hessian times the input vector.
-  success - flag indicating whether a solution was found.

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeHessian(), TaoLMVMSetSize()
@*/
int TaoAppHessianSolve(TAO_APPLICATION taoapp, Vec Vin, Vec Vout, PetscTruth *success){
  int     info;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  PetscValidHeaderSpecific(Vin,VEC_COOKIE,2);
  PetscValidHeaderSpecific(Vout,VEC_COOKIE,3);
  PetscStackPush("TAO Apply User Approximate Hessian");
  info = PetscLogEventBegin(Tao_HessianEval,taoapp,Vin,Vout,0);
  if (taoapp->hessiansolve){
    info = (*taoapp->hessiansolve)(taoapp,Vin,Vout,success,taoapp->usrhhhctx);
    CHKERRQ(info);
  }
  info = PetscLogEventEnd(Tao_HessianEval,taoapp,Vin,Vout,0);
  taoapp->nlsolve++;
  PetscStackPop;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppResetCounters"
/*@
   TaoAppResetCounters - Resent function evaluations counters to zero.

   Collective on TAO_APPLICATION

   Input Parameters:
.  taopp - the TAO_APPLICATION context

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeHessian
@*/
int TaoAppResetCounters(TAO_APPLICATION taoapp){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->nfeval=0;
  taoapp->ngeval=0;
  taoapp->nheval=0;
  taoapp->nlsolve=0;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppCounters"
/*@
   TaoAppCounters - Count the number of function, gradient, and
   Hessian evaluations, and the number of linear solves.

   Collective on TAO_APPLICATION

   Input Parameters:
.  taopp - the TAO_APPLICATION context

   Output Parameters:
   stat[4] ,the number of function, gradient, and Hessian evaluations. And
    the number of linear solves.

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeHessian
@*/
int TaoAppCounters(TAO_APPLICATION taoapp,TaoInt stats[4]){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  stats[0]=taoapp->nfeval;
  stats[1]=taoapp->ngeval;
  stats[2]=taoapp->nheval;
  stats[3]=taoapp->nlsolve;
  PetscFunctionReturn(0);
}


