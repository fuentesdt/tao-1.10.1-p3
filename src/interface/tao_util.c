/*$Id$*/

#include "src/tao_impl.h"       /*I   "tao_solver.h"   I*/


#undef __FUNCT__  
#define __FUNCT__ "TaoVecViewMonitor"
/*@C
   TaoVecViewMonitor - Monitors progress of the TAO_SOLVER solvers by calling 
   VecView() for the approximate solution at each iteration.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  dummy - either a viewer or TAO_NULL

   Level: advanced

.keywords: TAO_SOLVER, vector, monitor, view

.seealso: TaoSetMonitor(), TaoDefaultMonitor(), VecView()
@*/
int TaoVecViewMonitor(TAO_SOLVER tao,void *dummy)
{
  int    info;
  TaoVec *xx;

  TaoFunctionBegin;
  info = TaoGetSolution(tao,&xx);CHKERRQ(info);
  info = xx->View();CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoVecViewMonitorUpdate"
/*@C
   TaoVecViewMonitorUpdate - Monitors progress of the TAO_SOLVER solvers by calling 
   VecView() for the UPDATE to the solution at each iteration.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  dummy - either a viewer or TAO_NULL

   Level: advanced

.keywords: TAO_SOLVER, vector, monitor, view

.seealso: TaoSetMonitor(), TaoDefaultMonitor(), VecView()
@*/
int TaoVecViewMonitorUpdate(TAO_SOLVER tao,void *dummy)
{
  int    info;
  TaoVec *xx;

  TaoFunctionBegin;
  info = TaoGetStepDirectionVector(tao,&xx);CHKERRQ(info);
  if (xx){  info = xx->View();CHKERRQ(info); }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoDefaultMonitor"
/*@C
   TaoDefaultMonitor - Default routine for monitoring progress of the 
   TAO_SOLVER solvers (default).

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  dummy - unused context

   Notes:
   The routine prints the function value and gradient norm at each iteration.

   Level: advanced

.keywords: TAO_SOLVER, default, monitor, norm

.seealso: TaoSetMonitor(), TaoVecViewMonitor()
@*/
int TaoDefaultMonitor(TAO_SOLVER tao,void *dummy)
{
  int info;
  TaoInt its;
  double fct,gnorm;

  TaoFunctionBegin;
  its=tao->iter;
  fct=tao->fc;
  gnorm=tao->norm;
  info=TaoPrintInt(tao,"iter = %d,",its); CHKERRQ(info);
  info=TaoPrintDouble(tao," Function value: %12.10e,",fct); CHKERRQ(info);
  info=TaoPrintDouble(tao,"  Residual: %12.10e \n",gnorm);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoDefaultSMonitor"
/*
     Default (short) TAO_SOLVER Monitor, same as TaoDefaultMonitor() except
  it prints fewer digits of the residual as the residual gets smaller.
  This is because the later digits are meaningless and are often 
  different on different machines; by using this routine different 
  machines will usually generate the same output.
*/
int TaoDefaultSMonitor(TAO_SOLVER tao,void *dummy)
{
  int info;
  TaoInt its;
  double  fct,gnorm;

  TaoFunctionBegin;
  its=tao->iter;
  fct=tao->fc;
  gnorm=tao->norm;

  if (gnorm > 1.e-6) {
    info=TaoPrintInt(tao,"iter = %d,",its); CHKERRQ(info);
    info=TaoPrintDouble(tao," Function value %g,",fct); CHKERRQ(info);
    info=TaoPrintDouble(tao," Residual: %7.6f \n",gnorm);CHKERRQ(info);
  } else if (gnorm > 1.e-11) {
    info=TaoPrintInt(tao,"iter = %d,",its); CHKERRQ(info);
    info=TaoPrintDouble(tao," Function value %g,",fct); CHKERRQ(info);
    info=TaoPrintStatement(tao," Residual: < 1.0e-6 \n");
  } else {
    info=TaoPrintInt(tao,"iter = %d,",its); CHKERRQ(info);
    info=TaoPrintDouble(tao," Function value %g,",fct); CHKERRQ(info);
    info=TaoPrintStatement(tao," Residual: < 1.0e-11 \n");
  }
  TaoFunctionReturn(0);
}


/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoConverged_MaxIts"
/*@ 
   TaoConverged_MaxIts - Determines whether the solver has reached maximum number
   of iterations. 

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  dummy - unused dummy context

   Level: developer

.seealso: TaoSetTolerances(),TaoGetTerminationReason(),TaoSetTerminationReason()
@*/
int TaoConverged_MaxIts(TAO_SOLVER tao,void *dummy)
{
  int info;
  TaoInt its, maxits=tao->max_its;
  TaoTerminateReason reason=TAO_CONTINUE_ITERATING;

  TaoFunctionBegin;
  info = TaoGetSolutionStatus(tao,&its,0,0,0,0,0);CHKERRQ(info);
  if (its >= maxits) {
    info = PetscInfo2(tao,"TaoConverged_Default: Exceeded maximum number of iterations: %d > %d\n",its,maxits); CHKERRQ(info);
    reason = TAO_DIVERGED_MAXITS;
    info = TaoSetTerminationReason(tao,reason); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoConverged_Default"
/*@ 
   TaoConverged_Default - Determines whether the solver should continue iterating
   or terminate. 

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  dummy - unused dummy context

   Output Parameter:
.  reason - for terminating

   Notes:
   This routine checks the residual in the optimality conditions, the 
   relative residual in the optimity conditions, the number of function
   evaluations, and the function value to test convergence.  Some
   solvers may use different convergence routines.

   Level: developer

.seealso: TaoSetTolerances(),TaoGetTerminationReason(),TaoSetTerminationReason()
@*/
int TaoConverged_Default(TAO_SOLVER tao,void *dummy)
{
  int info;
  TaoInt its;
  TaoInt nfuncs=tao->nfuncs+tao->ngrads+tao->nfgrads, max_funcs=tao->max_funcs;
  double gnorm, gnorm0=tao->norm0;
  double f, trtol=tao->trtol,trradius=tao->step;
  double gatol,grtol,gttol,fatol,frtol,catol,crtol;
  double fmin=tao->fmin, cnorm, cnorm0=tao->cnorm0;
  double gnorm2;
  TaoTerminateReason reason=TAO_CONTINUE_ITERATING;

  TaoFunctionBegin;
  info = TaoGetSolutionStatus(tao,&its,&f,&gnorm,&cnorm,&trradius,&reason);
  info = TaoGetTolerances(tao,&fatol,&frtol,&catol,&crtol);CHKERRQ(info);
  info = TaoGetGradientTolerances(tao,&gatol,&grtol,&gttol);CHKERRQ(info);
  gnorm2=gnorm*gnorm;

  if (f != f ) {
    info = PetscInfo(tao,"TaoConverged_Default: Failed to converged, function is NaN\n"); CHKERRQ(info);
    reason = TAO_DIVERGED_NAN;
  } else if (f <= fmin && cnorm <=catol) {
    info = PetscInfo2(tao,"TaoConverged_Default: Converged due to function value %g < minimum function value %g\n", f,fmin); CHKERRQ(info);
    reason = TAO_CONVERGED_MINF;
  } else if (gnorm2 <= fatol && cnorm <=catol) {
    info = PetscInfo2(tao,"TaoConverged_Default: Converged due to residual norm %g < %g\n",gnorm2,fatol); CHKERRQ(info);
    reason = TAO_CONVERGED_ATOL;
  } else if (gnorm2 / TaoAbsScalar(f+1.0e-10)<= frtol && cnorm/TaoMax(cnorm0,1.0) <= crtol) {
    info = PetscInfo2(tao,"TaoConverged_Default: Converged due to relative residual norm %g < %g\n",gnorm2/TaoAbsScalar(f+1.0e-10),frtol); CHKERRQ(info);
    reason = TAO_CONVERGED_RTOL;
  } else if (gnorm<= gatol && cnorm <=catol) {
    info = PetscInfo2(tao,"TaoConverged_Default: Converged due to residual norm %g < %g\n",gnorm,gatol); CHKERRQ(info);
    reason = TAO_CONVERGED_ATOL;
  } else if ( f!=0 && TaoAbsScalar(gnorm/f) <= grtol && cnorm <= crtol) {
    info = PetscInfo3(tao,"TaoConverged_Default: Converged due to residual norm %g < |%g| %g\n",gnorm,f,grtol); CHKERRQ(info);
    reason = TAO_CONVERGED_ATOL;
  } else if (gnorm/gnorm0 <= gttol && cnorm <= crtol) {
    info = PetscInfo2(tao,"TaoConverged_Default: Converged due to relative residual norm %g < %g\n",gnorm/gnorm0,gttol); CHKERRQ(info);
    reason = TAO_CONVERGED_RTOL;
  } else if (nfuncs > max_funcs){
    info = PetscInfo2(tao,"TaoConverged_Default: Exceeded maximum number of function evaluations: %d > %d\n", nfuncs,max_funcs); CHKERRQ(info);
    reason = TAO_DIVERGED_MAXFCN;
  } else if ( tao->lsflag != 0 ){
    info = PetscInfo(tao,"TaoConverged_Default: Tao Line Search failure.\n"); CHKERRQ(info);
    reason = TAO_DIVERGED_LS_FAILURE;
  } else if (trradius < trtol && its > 0){
    info = PetscInfo2(tao,"TaoConverged_Default: Trust region/step size too small: %g < %g\n", trradius,trtol); CHKERRQ(info);
    reason = TAO_CONVERGED_TRTOL;
  } else {
    reason = TAO_CONTINUE_ITERATING;
  }
  info = TaoSetTerminationReason(tao,reason); CHKERRQ(info);
  info=TaoConverged_MaxIts(tao,0); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoCheckFG"
/*@C 
   TaoCheckFG - Check if the function and gradient vectors have
   been set properly and are compatible. 

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER context

   Level: developer

.seealso: TaoCheckFGH()
@*/
int TaoCheckFG(TAO_SOLVER tao)
{
  TaoVec *xx,*gg;
  TaoTruth flag;
  int info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);

  info=TaoGetSolution(tao,&xx);CHKERRQ(info);
  info=TaoGetGradient(tao,&gg);CHKERRQ(info);

  info=xx->Compatible(gg,&flag);CHKERRQ(info);
  if (flag==TAO_FALSE){
    SETERRQ(1,"Gradient and Variable vectors must have identical structure.");
  }

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCheckFGH"
/*@C 
   TaoCheckFGH - Check if the function, gradient vector, and
   Hessian matrix have been set properly and are compatible. 

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER context

   Level: developer

.seealso: TaoCheckFG()
@*/
int TaoCheckFGH(TAO_SOLVER tao)
{
  int info;
  TaoMat *HH;
  TaoFunctionBegin;
  
  info=TaoCheckFG(tao);CHKERRQ(info);
  info = TaoGetHessian(tao,&HH);CHKERRQ(info);
  if (!HH) {
    SETERRQ(1,"Must Provide Hessian Matrix");
  }

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoCheckConstraints"
/*@C 
   TaoCheckConstraints - Check if the nonlinear constraints have
   been set and if the data structures are compatible.

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER context

   Level: developer

.seealso: TaoCheckFG()
@*/
int TaoCheckConstraints(TAO_SOLVER tao)
{
  TaoVec *solu, *cons;
  TaoMat *jac;
  TaoTruth flag;
  int info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);

  info = TaoGetSolution(tao, &solu); CHKERRQ(info);
  info = TaoGetConstraints(tao, &cons); CHKERRQ(info);
  info = TaoGetJacobian(tao, &jac); CHKERRQ(info);
  
  info = jac->Compatible(solu,cons,&flag);CHKERRQ(info);
  if (flag == TAO_FALSE){
    SETERRQ(1,"Jacobian matrix not consistent with Variable Vector or Constraint Vector");
  }

  TaoFunctionReturn(0);
}
  
#undef __FUNCT__  
#define __FUNCT__ "TaoCheckBounds"
/*@C 
   TaoCheckBounds - Check if the variable bounds have been
   set and the data structures are compatible.

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER context

   Level: developer

.seealso: TaoCheckFG()
@*/
int TaoCheckBounds(TAO_SOLVER tao)
{
  int info;
  TaoTruth flag;
  TaoVec *xx,*xxll,*xxuu;

  info=TaoGetSolution(tao,&xx);CHKERRQ(info);
  info=TaoGetVariableBounds(tao,&xxll,&xxuu);CHKERRQ(info);

  info=xx->Compatible(xxll,&flag);CHKERRQ(info);
  if (flag==TAO_FALSE){
    SETERRQ(1,"Vector of lower bounds not Compatible with Variable Vector");
  }

  info=xx->Compatible(xxuu,&flag);CHKERRQ(info);
  if (flag==TAO_FALSE){
    SETERRQ(1,"Vector of upper bounds not Compatible with Variable vector");
  }

  TaoFunctionReturn(0);
}

static int TaoDefaultMeritFunction(TAO_SOLVER,TaoVec*,double*,void*);
static int TaoDefaultMeritFunctionGradient(TAO_SOLVER,TaoVec*,double*,TaoVec*,void*);
static int TaoDefaultMeritGradient(TAO_SOLVER,TaoVec*,TaoVec*,void*);

#undef __FUNCT__  
#define __FUNCT__ "TaoDefaultMeritFunction"
static int TaoDefaultMeritFunction(TAO_SOLVER tao, TaoVec *xx, double *f, void *ctx)
{
  int info;

  TaoFunctionBegin;
  info=TaoComputeFunction(tao,xx,f);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoDefaultMeritFunctionGradient"
static int TaoDefaultMeritFunctionGradient(TAO_SOLVER tao, TaoVec *xx, double *f, TaoVec *gg, void *ctx)
{
  int info;

  TaoFunctionBegin;
  info=TaoComputeFunctionGradient(tao,xx,f,gg);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoDefaultMeritGradient"
static int TaoDefaultMeritGradient(TAO_SOLVER tao, TaoVec *xx, TaoVec *gg, void *ctx)
{
  int info;

  TaoFunctionBegin;
  info=TaoComputeGradient(tao,xx,gg);CHKERRQ(info);
  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoSetDefaultMeritFunction"
/*@
   TaoSetDefaultMeritFunction - Set the merit function equal to the
   objective function

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER context

   Level: developer

@*/
int TaoSetDefaultMeritFunction(TAO_SOLVER tao)
{
  int info;
  TaoFunctionBegin;
  info=TaoSetMeritFunction(tao,TaoDefaultMeritFunction,TaoDefaultMeritFunctionGradient,TaoDefaultMeritGradient,0,0,0);
  CHKERRQ(info);
  TaoFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "TaoSetQuadraticMeritFunction"
/* @C 
   TaoSetQuadraticMeritFunction - Set the merit function to a quadratic
   objective function 1/2 x^T A x + b^T x + c

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
.  A - The matrix term of the quadratic function
.  b - The linear term of the quadratic function
-  c - The constant term of the quadratic funtion

   Level: developer

@ */
int TaoSetQuadraticMeritFunction(TAO_SOLVER tao, TaoMat *A, TaoVec*b, double c)
{
  TaoFunctionBegin;
  SETERRQ(1,"TAO ERROR: Not implemented");
  //  TaoFunctionReturn(0);
}









