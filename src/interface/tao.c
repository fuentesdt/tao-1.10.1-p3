/*$Id$*/

#include "src/tao_impl.h"      /*I "tao_solver.h"  I*/


#undef __FUNCT__  
#define __FUNCT__ "TaoMonitor"
/*@C
  TaoMonitor - Monitor the solver and the current solution.  This
  routine will calls the records the iteration number and residual statistics,
  monitors specified by the user, and calls the termaination routine.

   Input Parameters:
+  tao - the TAO_SOLVER context
.  f - the current objective function value
.  iterate - the current iterate number (>=0)
.  fnorm - the gradient norm, square root of the duality gap, or other measure
indicating distince from optimality.  This measure will be recorded and
used for some termination tests.
.  cnorm - the infeasibility of the current solution with regard to the constraints.
-  step - multiple of the step direction added to the previous iterate.

   Output Parameters:
.  reason - The termination reason, which can equal TAO_CONTINUE_ITERATING

   Options Database Key:
.  -tao_monitor - The default monitor, which prints statistics to standard output is used.

.seealso TaoGetTerminationReason(),TaoGetSolutionStatus()

   Level: developer

.keywords: Monitor, convergence
@*/
int TaoMonitor(TAO_SOLVER tao, TaoInt iterate, double f, double fnorm, double cnorm, double step, TaoTerminateReason *reason)
{
  int info;
  TaoInt        i;
  TaoTruth cstop;
  TaoFunctionBegin;
  if (iterate>=tao->iter){
    TaoLogConvHistory(tao,fnorm,tao->iter);
    tao->norm=fnorm; tao->cnorm=cnorm; tao->fc=f;
    tao->iter=TaoMax(tao->iter,iterate);
    tao->step=step;
  }
  if (iterate==0){ 
    tao->norm0=fnorm; tao->cnorm0=cnorm;
  }
  info=TaoCheckConvergence(tao,&tao->reason);CHKERRQ(info);
  if (iterate>0){
    info=tao->taoappl->Monitor2(tao->vec_sol,tao->vec_grad,tao->vec_sol_update,&cstop);CHKERRQ(info);
    if (cstop==TAO_TRUE) tao->reason=TAO_CONVERGED_USER;
  }
  for ( i=0; i<tao->numbermonitors; i++ ) {
    info = (*tao->monitor[i])(tao,tao->monitorcontext[i]);CHKERRQ(info);
  }
  info=tao->taoappl->Monitor();CHKERRQ(info);
  *reason = tao->reason;

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetSolutionStatus"
/*@C
  TaoGetSolutionStatus - Get the current iterate, objective value, residual, 
  infeasibility, and termination 

   Input Parameters:
.  tao - the TAO_SOLVER context

   Output Parameters:
+  iterate - the current iterate number (>=0)
.  f - the current function value
.  gnorm - the square of the gradient norm, duality gap, or other measure
indicating distance from optimality.
.  cnorm - the infeasibility of the current solution with regard to the constraints.
.  xdiff - the step length or trust region radius of the most recent iterate.
-  reason - The termination reason, which can equal TAO_CONTINUE_ITERATING

   Level: intermediate

   Note:
   TAO returns the values set by the solvers in the routine TaoMonitor().

   Note:
   If any of the output arguments are set to TAO_NULL, no value will be 
   returned.


.seealso: TaoMonitor(), TaoGetTerminationReason()

.keywords: convergence, monitor
@*/
int TaoGetSolutionStatus(TAO_SOLVER tao, TaoInt* iterate, double* f, double* gnorm, double *cnorm, double *xdiff, TaoTerminateReason *reason)
{

  TaoFunctionBegin;
  if (iterate) *iterate=tao->iter;
  if (f) *f=tao->fc;
  if (gnorm) *gnorm=tao->norm;
  if (cnorm) *cnorm=tao->cnorm;
  if (reason) *reason=tao->reason;
  if (xdiff) *xdiff=tao->step;

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoCheckConvergence"
/*@C

   TaoCheckConvergence - Checks the convergence of the solver

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER context

   Output Parameters:
.  reason - one of

$    TAO_CONVERGED_ATOL (2),        (res <= atol)  
$    TAO_CONVERGED_RTOL (3),        (res/res0 <= rtol) 
$    TAO_CONVERGED_TRTOL (4),       (xdiff <= trtol) 
$    TAO_CONVERGED_MINF (5),        (f <= fmin)
$    TAO_CONVERGED_USER (6),        (user defined)

$    TAO_DIVERGED_MAXITS (-2),      (its>maxits)
$    TAO_DIVERGED_NAN (-4),         (Numerical problems)
$    TAO_DIVERGED_MAXFCN (-5),      (nfunc > maxnfuncts)
$    TAO_DIVERGED_LS_FAILURE (-6),  (line search failure)
$    TAO_DIVERGED_TR_REDUCTION (-7),
$    TAO_DIVERGED_USER (-8),        (user defined)

$    TAO_CONTINUE_ITERATING  (0)

   where
+  res - residual of optimality conditions
.  res0 - initial residual of optimality conditions
.  xdiff - current trust region size
.  f - function value
.  atol - absolute tolerance
.  rtol - relative tolerance
.  its - current iterate number
.  maxits - maximum number of iterates
.  nfunc - number of function evaluations
-  maxnfuncts - maximum number of function evaluations


   Level: advanced

.seealso: TaoGetTerminationReason()

.keywords: Convergence

@*/
int TaoCheckConvergence(TAO_SOLVER tao, TaoTerminateReason *reason)
{
  int        info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (tao->converged){
    info = (*tao->converged)(tao,tao->cnvP);CHKERRQ(info);
  }
  if (reason) *reason=tao->reason;

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoView"
/*@ 
   TaoView - Prints the TAO_SOLVER data structure.

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER context

   Options Database Key:
.  -tao_view - Calls TaoView() at end of TaoSolve()

   Level: beginner

.seealso: TaoGetTerminationReason(), TaoGetSolutionStatus(), TaoViewLinearSolver()

.keywords: View

@*/
int TaoView(TAO_SOLVER tao)
{
  int        info;
  TaoMethod  type;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);

  info = TaoPrintStatement(tao,"TAO_SOLVER:\n");CHKERRQ(info);
  info = TaoGetMethod(tao,&type);CHKERRQ(info);
  if (type) {
    info = TaoPrintString(tao,"  method: %s\n",type);CHKERRQ(info);
  } else {
    info = TaoPrintStatement(tao,"  method: not set yet\n");CHKERRQ(info);
  }
  if (tao->view) {
    info = (*tao->view)(tao,tao->data);CHKERRQ(info);
  }
  
  info=TaoPrintDouble(tao,"  convergence tolerances: fatol=%g,",tao->fatol);CHKERRQ(info);
  info=TaoPrintDouble(tao," frtol=%g\n",tao->frtol);CHKERRQ(info);

  info=TaoPrintDouble(tao,"  convergence tolerances: gatol=%g,",tao->gatol);CHKERRQ(info);
  info=TaoPrintDouble(tao," trtol=%g,",tao->trtol);CHKERRQ(info);
  info=TaoPrintDouble(tao," gttol=%g\n",tao->gttol);CHKERRQ(info);

  info = TaoPrintDouble(tao,"  Residual in Function/Gradient:=%e\n",tao->norm);CHKERRQ(info);

  if (tao->cnorm>0 || tao->catol>0 || tao->crtol>0){
    info=TaoPrintStatement(tao,"  convergence tolerances:");CHKERRQ(info);
    info=TaoPrintDouble(tao," catol=%g,",tao->catol);CHKERRQ(info);
    info=TaoPrintDouble(tao," crtol=%g\n",tao->crtol);CHKERRQ(info);
    info = TaoPrintDouble(tao,"  Residual in Constraints:=%e\n",tao->cnorm);CHKERRQ(info);
  }

  if (tao->trtol>0){
    info=TaoPrintDouble(tao,"  convergence tolerances: trtol=%g\n",tao->trtol);CHKERRQ(info);
    info=TaoPrintDouble(tao,"  Final step size/trust region radius:=%g\n",tao->step);CHKERRQ(info);
  }

  if (tao->fmin>-1.e25){
    info=TaoPrintDouble(tao,"  convergence tolerances: function minimum=%g\n",tao->fmin);CHKERRQ(info);
  }
  info = TaoPrintDouble(tao,"  Objective value=%e\n",tao->fc);CHKERRQ(info);

  info = TaoPrintInt(tao,"  total number of iterations=%d,          ",tao->iter);CHKERRQ(info);
  info = TaoPrintInt(tao,"              (max: %d)\n",tao->max_its);CHKERRQ(info);

  if (tao->nfuncs>0){
    info = TaoPrintInt(tao,"  total number of function evaluations=%d,",tao->nfuncs);CHKERRQ(info);
    info = TaoPrintInt(tao,"                max: %d\n",tao->max_funcs);CHKERRQ(info);
  }
  if (tao->ngrads>0){
    info = TaoPrintInt(tao,"  total number of gradient evaluations=%d,",tao->ngrads);CHKERRQ(info);
    info = TaoPrintInt(tao,"                max: %d\n",tao->max_funcs);CHKERRQ(info);
  }
  if (tao->nfgrads>0){
    info = TaoPrintInt(tao,"  total number of function/gradient evaluations=%d,",tao->nfgrads);CHKERRQ(info);
    info = TaoPrintInt(tao,"    (max: %d)\n",tao->max_funcs);CHKERRQ(info);
  }
  if (tao->nhesss>0){
    info = TaoPrintInt(tao,"  total number of Hessian evaluations=%d\n",tao->nhesss);CHKERRQ(info);
  }
  if (tao->linear_its>0){
    info = TaoPrintInt(tao,"  total Krylov method iterations=%d\n",tao->linear_its);CHKERRQ(info);
  }
  if (tao->nvfunc>0){
    info = TaoPrintInt(tao,"  total number of constraint function evaluations=%d\n",tao->nvfunc);CHKERRQ(info);
  }
  if (tao->njac>0){
    info = TaoPrintInt(tao,"  total number of Jacobian evaluations=%d\n",tao->njac);CHKERRQ(info);
  }

  if (tao->reason>0){
    info = TaoPrintStatement(tao,"  Solution found\n");CHKERRQ(info);
  } else {
    info = TaoPrintInt(tao,"  Solver terminated: %d\n",tao->reason);CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}



/* ----------- Routines to set solver parameters ---------- */

#undef __FUNCT__  
#define __FUNCT__ "TaoSetGradientTolerances"
/*@
   TaoSetGradientTolerances - Sets the stopping criteria in terms of the norm
   of the Lagrangian function.  The algorithm will terminate when the norm
   of the gradient is less that the absolute tolerance, or when the norm
   of the gradient has been reduced by a factor of the reduction tolerance, 
   or when the norm of the gradient divided by the absolute value of the 
   objective function is less than the relative tolerance.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
.  gatol - stop if norm of gradient is less than 
.  grtol - stop if relative norm of gradient is less than
-  gttol - stop if norm of gradient is reduced by a factor of

   Options Database Keys: 
+  -tao_gatol <gatol> - sets gatol
.  -tao_grtol <grtol> - sets grtol
-  -tao_gttol <gttol> - sets gttol

   Level: intermediate

.keywords: Gradient, options, convergence

.seealso: TaoSetTolerances(), TaoGetGradientTolerances()
@*/
int TaoSetGradientTolerances(TAO_SOLVER tao,double gatol, double grtol, double gttol)
{
  int info;
  double zero=0.0;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);

  if (gatol==TAO_DEFAULT){}
  else if (gatol<0){
    info=PetscInfo(tao,"Absolute Gradient tolerance < 0 and ignored"); CHKERRQ(info);
    CHKERRQ(info); } 
  else{
    tao->gatol      = TaoMax(zero,gatol);
  } 

  if (grtol==TAO_DEFAULT){}
  else if (grtol<0){
    info=PetscInfo(tao,"Relative Gradient tolerance < 0 and ignored");
    CHKERRQ(info); } 
  else{
    tao->grtol      = TaoMax(zero,grtol);
  } 

  if (gttol==TAO_DEFAULT){}
  else if (gttol<0){
    info=PetscInfo(tao,"Gradient reduction tolerance < 0 and ignored");
    CHKERRQ(info); } 
  else{
    tao->gttol      = TaoMax(zero,gttol);
  } 


  TaoFunctionReturn(0);
}

/* ----------- Routines to set solver parameters ---------- */

#undef __FUNCT__  
#define __FUNCT__ "TaoGetGradientTolerances"
/*@
   TaoGetGradientTolerances - Returns the gradient termination tolerances.

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER tao context

   Output Parameters:
+  gatol - the absolute gradient tolerance
.  grtol - the relative gradient tolerance
-  gttol - the gradient reduction tolerance

   Level: intermediate

.keywords: options, convergence, View

.seealso: TaoSetGradientTolerances()
@*/
int TaoGetGradientTolerances(TAO_SOLVER tao,double *gatol, double *grtol, double *gttol)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (gatol) *gatol = tao->gatol;
  if (grtol) *grtol = tao->grtol;
  if (gttol) *gttol = tao->gttol;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetFunctionLowerBound"
/*@
   TaoSetFunctionLowerBound - Sets a bound on the solution objective value.
   When an approximate solution with an objective value below this number
   has been found, the solver will terminate.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
-  fmin - the tolerance

   Options Database Keys: 
.    -tao_fmin <fmin> - sets the minimum function value

   Level: intermediate

.keywords: options, View, Bounds,

.seealso: TaoSetTolerances()
@*/
int TaoSetFunctionLowerBound(TAO_SOLVER tao,double fmin)
{
  double dflt=-1.0e+30;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (fmin != TAO_DEFAULT)  tao->fmin = fmin;
  else tao->fmin=dflt;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetMaximumIterates"
/*@
   TaoSetMaximumIterates - Sets a maximum number of iterates.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
-  maxits - the maximum number of iterates (>=0)

   Options Database Keys: 
.    -tao_max_its <its> - sets the maximum number of iterations

   Level: intermediate

.keywords: options, Iterate, convergence

.seealso: TaoSetTolerances(), TaoSetMaximumFunctionEvaluations()
@*/
int TaoSetMaximumIterates(TAO_SOLVER tao,TaoInt maxits)
{
  TaoInt zero=0;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (maxits != TAO_DEFAULT)  tao->max_its = TaoMax(zero,maxits);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetMaximumFunctionEvaluations"
/*@
   TaoSetMaximumFunctionEvaluations - Sets a maximum number of 
   function evaluations.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
-  nfcn - the maximum number of function evaluations (>=0)

   Options Database Keys: 
.    -tao_max_funcs <nfcn> - sets the maximum number of function evaluations

   Level: intermediate

.keywords: options, Iterate,  convergence

.seealso: TaoSetTolerances(), TaoSetMaximumIterates()
@*/
int TaoSetMaximumFunctionEvaluations(TAO_SOLVER tao,TaoInt nfcn)
{
  int zero=0;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (nfcn != TAO_DEFAULT)  tao->max_funcs = TaoMax(zero,nfcn);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetTolerances"
/*@
   TaoSetTolerances - Sets convergence parameters. TAO tries to satisfy an
   absolute stopping criteria or a relative stopping criteria.
   
   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
.  fatol - absolute convergence tolerance
.  frtol - relative convergence tolerance
.  catol -  allowable error in constraints
-  crtol - allowable relative error in constraints

   Options Database Keys: 
+  -tao_fatol <fatol> - Sets fatol
.  -tao_frtol <frtol> - Sets frtol
.  -tao_catol <catol> - Sets catol
-  -tao_crtol <crtol> - Sets crtol

   Absolute Stopping Criteria:
$  f <= f + fatol
$  B1 - catol <= B(X) <= B2 + catol

   Relative stopping criteria:
$  f <= f + frtol*|f|
$  B1 - catol <= B(X) <= B2 + catol

   Level: beginner

.keywords: options, convergence

.seealso: TaoSetMaximumIterates(),TaoSetTrustRegionTolerance(), TaoSetGradientTolerances
TaoSetMaximumFunctionEvaluations()
@*/
int TaoSetTolerances(TAO_SOLVER tao,double fatol,double frtol,double catol,double crtol)
{
  int info;
  double zero=0.0;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);

  if (fatol==TAO_DEFAULT){}
  else if (fatol<0){
    info=PetscInfo(tao,"Absolute convergence tolerance < 0 and ignored");
    CHKERRQ(info); } 
  else{
    tao->fatol      = TaoMax(zero,fatol);
  } 

  if (frtol==TAO_DEFAULT){}
  else if (frtol<0){
    info=PetscInfo(tao,"Relative convergence tolerance < 0 and ignored");
    CHKERRQ(info); } 
  else{
    tao->frtol      = TaoMax(zero,frtol);
  } 

  if (catol==TAO_DEFAULT){}
  else if (catol<0){
    info=PetscInfo(tao,"Absolute constraint tolerance < 0 and ignored");
    CHKERRQ(info); } 
  else{
    tao->catol      = TaoMax(zero,catol);
  } 

  if (crtol==TAO_DEFAULT){}
  else if (crtol<0){
    info=PetscInfo(tao,"Relative constraint tolerance < 0 and ignored");
    CHKERRQ(info); } 
  else{
    tao->crtol      = TaoMax(zero,crtol);
  } 

  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoGetTolerances"
/*@
   TaoGetTolerances - Gets convergence parameters.
   

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER solver context

   Input Parameters:
+  fatol - absolute convergence tolerance
.  frtol - relative convergence tolerance
.  catol -  trust region convergence tolerance
-  crtol - convergence of the function evaluates less than this tolerance

   Level: advanced

.keywords: options, convergence

.seealso: TaoSetTolerances
@*/
int TaoGetTolerances(TAO_SOLVER tao,double *fatol,double *frtol,double *catol,double *crtol)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (fatol) *fatol = tao->fatol;
  if (frtol) *frtol = tao->frtol;
  if (catol) *catol = tao->catol;
  if (crtol) *crtol = tao->crtol;
  TaoFunctionReturn(0);
}


/* ------------ Routines to set performance monitoring options ----------- */

#undef __FUNCT__  
#define __FUNCT__ "TaoSetMonitor"
/*@C
   TaoSetMonitor - Sets an ADDITIONAL function that is to be used at every
   iteration of the unconstrained minimization solver to display the iteration's 
   progress.   

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
.  mymonitor - monitoring routine
-  mctx - [optional] user-defined context for private data for the 
          monitor routine (may be TAO_NULL)

   Calling sequence of mymonitor:
$     int mymonitor(TAO_SOLVER tao,void *mctx)

+    tao - the TAO_SOLVER solver context
-    mctx - [optional] monitoring context


   Options Database Keys:
+    -tao_monitor        - sets TaoDefaultMonitor()
.    -tao_smonitor        - sets short monitor
-    -tao_cancelmonitors - cancels all monitors that have been hardwired into a code by calls to TaoSetMonitor(), but does not cancel those set via the options database.

   Notes: 
   Several different monitoring routines may be set by calling
   TaoSetMonitor() multiple times; all will be called in the 
   order in which they were set.

   Level: intermediate

.keywords: options, monitor, View

.seealso: TaoDefaultMonitor(), TaoClearMonitor(),  TaoSetDestroyRoutine()
@*/
int TaoSetMonitor(TAO_SOLVER tao,int (*mymonitor)(TAO_SOLVER,void*),void *mctx)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (tao->numbermonitors >= MAX_TAO_MONITORS) {
    SETERRQ(1,"Too many monitors set");
  }

  tao->monitor[tao->numbermonitors]           = mymonitor;
  tao->monitorcontext[tao->numbermonitors++]  = (void*)mctx;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoClearMonitor"
/*@
   TaoClearMonitor - Clears all the monitor functions for a TAO_SOLVER object.

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER solver context

   Options Database:
.  -tao_cancelmonitors - cancels all monitors that have been hardwired
    into a code by calls to TaoSetMonitor(), but does not cancel those 
    set via the options database

   Notes: 
   There is no way to clear one specific monitor from a TAO_SOLVER object.

   Level: advanced

.keywords: options, monitor, View

.seealso: TaoDefaultMonitor(), TaoSetMonitor()
@*/
int TaoClearMonitor(TAO_SOLVER tao)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  tao->numbermonitors = 0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetConvergenceTest"
/*@C
   TaoSetConvergenceTest - Sets the function that is to be used 
   to test for convergence of the iterative minimization solution.  The 
   new convergence testing routine will replace TAO's default 
   convergence test.   

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
.  conv - routine to test for convergence
-  cctx - [optional] context for private data for the convergence routine 
          (may be TAO_NULL)

   Calling sequence of conv:
$     int conv (TAO_SOLVER tao, void *cctx)

+    tao - the TAO_SOLVER solver context
-    cctx - [optional] convergence context

   Note: The new convergence testing routine should call TaoSetTerminationReason().

   Level: intermediate

.keywords: options, convergence

.seealso: TaoSetTerminationReason(), TaoGetSolutionStatus(), TaoGetTolerances(),  TaoGetGradientTolerances()

@*/
int TaoSetConvergenceTest(TAO_SOLVER tao,int (*conv)(TAO_SOLVER,void*),void *cctx)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  (tao)->converged = conv;
  (tao)->cnvP      = cctx;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetTerminationReason"
/*@C
   TaoGetTerminationReason - Gets the reason the TAO_SOLVER iteration was stopped.

   Not Collective

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Output Parameter:
.  reason - one of

$    TAO_CONVERGED_ATOL (2),        (res <= atol)  
$    TAO_CONVERGED_RTOL (3),        (res/res0 <= rtol) 
$    TAO_CONVERGED_TRTOL (4),       (xdiff <= trtol) 
$    TAO_CONVERGED_MINF (5),        (f <= fmin)
$    TAO_CONVERGED_USER (6),        (user defined)

$    TAO_DIVERGED_MAXITS (-2),      (its>maxits)
$    TAO_DIVERGED_NAN (-4),         (Numerical problems)
$    TAO_DIVERGED_MAXFCN (-5),      (nfunc > maxnfuncts)
$    TAO_DIVERGED_LS_FAILURE (-6),  (line search failure)
$    TAO_DIVERGED_TR_REDUCTION (-7),
$    TAO_DIVERGED_USER (-8),        (user defined)

$    TAO_CONTINUE_ITERATING  (0)

   where
+  res - residual of optimality conditions
.  res0 - initial residual of optimality conditions
.  xdiff - current trust region size
.  f - function value
.  atol - absolute tolerance
.  rtol - relative tolerance
.  its - current iterate number
.  maxits - maximum number of iterates
.  nfunc - number of function evaluations
-  maxnfuncts - maximum number of function evaluations

   Level: intermediate

.keywords: convergence, View

.seealso: TaoSetConvergenceTest(), TaoSetTolerances()
@*/
int TaoGetTerminationReason(TAO_SOLVER tao,TaoTerminateReason *reason)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  *reason = tao->reason;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetConvergenceHistory"
/*@
   TaoSetConvergenceHistory - Sets the array used to hold the convergence history.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
.  a   - array to hold history
.  its - integer array holds the number of linear iterations (or
         negative if not converged) for each solve.
.  na  - size of a and its
-  reset - TAO_TRUE indicates each new minimization resets the history counter to zero,
           else it continues storing new values for new minimizations after the old ones

   Notes:
   If set, this array will contain the gradient norms computed at each step.

   This routine is useful, e.g., when running a code for purposes
   of accurate performance monitoring, when no I/O should be done
   during the section of code that is being timed.

   Level: intermediate

.keywords: options, view, monitor, convergence, history

.seealso: TaoGetConvergenceHistory()

@*/
int TaoSetConvergenceHistory(TAO_SOLVER tao, double *a, TaoInt *its,TaoInt na,TaoTruth reset)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (na) TaoValidScalarPointer(a,2);
  tao->conv_hist       = a;
  tao->conv_hist_its   = its;
  tao->conv_hist_max   = na;
  tao->conv_hist_reset = reset;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetConvergenceHistory"
/*@C
   TaoGetConvergenceHistory - Gets the array used to hold the convergence history.

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Output Parameters:
+  a   - array to hold history
.  its - integer array holds the number of linear iterations (or
         negative if not converged) for each solve.
-  na  - size of a and its

   Notes:
    The calling sequence for this routine in Fortran is
$   call TaoGetConvergenceHistory(TAO_SOLVER tao, integer na, integer info)

   This routine is useful, e.g., when running a code for purposes
   of accurate performance monitoring, when no I/O should be done
   during the section of code that is being timed.

   Level: advanced

.keywords: convergence, history, monitor, View

.seealso: TaoSetConvergencHistory()

@*/
int TaoGetConvergenceHistory(TAO_SOLVER tao, double **a, TaoInt **its,TaoInt *na)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (a)   *a   = tao->conv_hist;
  if (its) *its = tao->conv_hist_its;
  if (na) *na   = tao->conv_hist_len;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSolve"
/*@
   TaoSolve - Solves an unconstrained minimization problem.  Call TaoSolve() 
   after calling TaoCreate() and optional routines of the form TaoSetXXX().

   Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER solver context

   Notes:
   By default the TAO solvers use an initial starting guess of zero.  To
   provide an alternative initial guess, the user must call TaoAppSetInitialSolutionVec()
   before calling TaoSolve().

   Level: advanced

.keywords: Solve

.seealso: TaoCreate(), TaoDestroy()
@*/
int TaoSolve(TAO_SOLVER tao)
{
  int      info;
  TaoVec   *xx;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = TaoGetSolution(tao,&xx);CHKERRQ(info);
  info = tao->taoappl->InitializeVariables(xx);CHKERRQ(info);
  info = TaoSetUp(tao);CHKERRQ(info);
  info = TaoSetDefaultStatistics(tao); CHKERRQ(info);
  if (tao->solve){ info = (*(tao)->solve)(tao,tao->data);CHKERRQ(info); }
  if (tao->viewtao) { info = TaoView(tao);CHKERRQ(info); }
  if (tao->viewksptao) { info = TaoViewLinearSolver(tao);CHKERRQ(info); }
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetMethod"
/*@C
   TaoGetMethod - Gets the TAO_SOLVER method type and name (as a string).

   Not Collective

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Output Parameter:
.  type - TAO_SOLVER method (a charactor string)

   Level: intermediate

.keywords: method
@*/
int TaoGetMethod(TAO_SOLVER tao, TaoMethod *type)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  *type = ((PetscObject)tao)->type_name;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetStepDirectionVector"
/*@C
   TaoSetStepDirectionVector - Sets the vector where the solution update is
   stored. 

   Not Collective, but Vec is parallel if TAO_SOLVER is parallel

   Input Parameter:
+  tao - the TAO_SOLVER solver context
-  dx - the step direction vector

   Level: developer

.keywords: solution, step direction

.seealso: TaoGetSolution(), TaoGetStepDirectionVector()
@*/
int TaoSetStepDirectionVector(TAO_SOLVER tao,TaoVec* dx)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  tao->vec_sol_update=dx;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetStepDirectionVector"
/*@C
   TaoGetStepDirectionVector - Returns the vector where the solution update is
   stored. 

   Not Collective, but Vec is parallel if TAO_SOLVER is parallel

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Output Parameter:
.  xx - the solution update

   Level: developer

.keywords: solution

.seealso: TaoGetSolution()
@*/
int TaoGetStepDirectionVector(TAO_SOLVER tao,TaoVec** xx)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  *xx = tao->vec_sol_update;
  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoSetTrustRegionTolerance"
/*@
  TaoSetTrustRegionTolerance - Sets a minimum step size or trust region radius.  The
  solver will terminate when the step size or radius of the trust region is smaller
  than this tolerance.
  
  Collective on TAO_SOLVER
  
  Input Parameters:
+  tao - the TAO_SOLVER solver context
-  steptol - tolerance
  
  Options Database Key: 
.  -tao_steptol <trtol> - Sets steptol
  
  Level: intermediate
  
.keywords: options, convergence, trust region
  
.seealso: TaoSetTolerances()
@*/
int TaoSetTrustRegionTolerance(TAO_SOLVER tao,double steptol)
{  
  double zero=0.0, dflt=0.0;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (steptol != TAO_DEFAULT)  tao->xtol = TaoMax(zero,steptol);
  else tao->xtol=dflt;
  tao->trtol=steptol;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetTrustRegionRadius"
/*@
   TaoGetTrustRegionRadius - Gets the current trust region radius

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - a TAO optimization solver

   Output Parameter:
.  radius - the trust region radius

   Level: advanced

.keywords: options, view, trust region

.seealso: TaoSetTrustRegionRadius(), TaoSetTrustRegionTolerance()
@*/
int TaoGetTrustRegionRadius(TAO_SOLVER tao,double *radius)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  *radius=tao->step;
  TaoFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "TaoGetInitialTrustRegionRadius"
/*@C
   TaoGetInitialTrustRegionRadius - Gets the initial trust region radius

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - a TAO optimization solver

   Output Parameter:
.  radius - the initial trust region radius

   Level: intermediate

.keywords: options, trust region

.seealso: TaoGetTrustRegionRadius(), TaoSetTrustRegionRadius()
@*/
int TaoGetInitialTrustRegionRadius(TAO_SOLVER tao,double *radius)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  *radius=tao->trust0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetTrustRegionRadius"
/*@
   TaoSetTrustRegionRadius - Sets the initial trust region radius.

   Collective on TAO_SOLVER

   Input Parameter:
+  tao - a TAO optimization solver
-  radius - the trust region radius

   Level: intermediate

   Options Database Key:
.  -tao_trust0

.keywords: trust region

.seealso: TaoGetTrustRegionRadius(), TaoSetTrustRegionTolerance()
@*/
int TaoSetTrustRegionRadius(TAO_SOLVER tao,double radius)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  tao->trust0=radius;
  tao->step=radius;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetVariableBounds"
/*@
   TaoSetVariableBounds - Sets lower and upper bounds on the variables.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
.  xxll - vector of lower bounds upon the solution vector
-  xxuu - vector of upper bounds upon the solution vector

   Level: developer

.keywords: bounds

.seealso: TaoGetVariableBounds(), TaoAppSetVariableBounds()
@*/
int TaoSetVariableBounds(TAO_SOLVER tao,TaoVec *xxll,TaoVec *xxuu)
{
  int info;
  double dd;
  TaoVec *xx;
  TaoTruth flag;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  info=TaoGetSolution(tao,&xx);CHKERRQ(info);
  if (xxll){
    info=xx->Compatible(xxll,&flag); CHKERRQ(info);
    if (flag == TAO_FALSE){
      SETERRQ(1,"Vector of lower bounds not Compatible with Variable Vector");
    }
    if (tao->XL==0){
      dd=-TAO_INFINITY;
      info = xxll->SetToConstant(dd); CHKERRQ(info);
    }
  }
  if (xxuu){
    info=xx->Compatible(xxuu,&flag); CHKERRQ(info);
    if (flag == TAO_FALSE){
      SETERRQ(1,"Vector of upper bounds not Compatible with Variable vector");
    }
    if (tao->XU==0){
      dd= TAO_INFINITY;
      info = xxuu->SetToConstant(dd); CHKERRQ(info);
    }
  } 
  tao->XL=xxll;
  tao->XU=xxuu;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetDualVariables"
/*@
   TaoGetDualVariables - Gets the dual variables corresponding to
   the bounds of the variables.

   Collective on TAO_SOLVER

   Input Parameter:
+  tao - the TAO_SOLVER solver context
.  DXL   - vector to place the dual variables of the lower bounds
-  DXU   - vector to place the dual variables of the upper bounds

   Output Parameter:
+  DXL   - dual variables of the lower bounds
-  DXU   - dual variables of the upper bounds

   Level: advanced

.keywords: dual, bounds

.seealso: TaoGetVariableBounds()
@*/
int TaoGetDualVariables(TAO_SOLVER tao, TaoVec *DXL, TaoVec *DXU)
{
  int info;
  TaoTruth flag;
  TaoVec *XL, *XU;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  info=TaoGetVariableBounds(tao,&XL,&XU);CHKERRQ(info);

  info=DXU->Compatible(XU,&flag); CHKERRQ(info);
  if (flag == TAO_FALSE){
    SETERRQ(1,"Dual bound vectors not Compatible");
  }

  info=DXL->Compatible(XL,&flag); CHKERRQ(info);
  if (flag == TAO_FALSE){
    SETERRQ(1,"Dual bound vectors not Compatible");
  }

  if (tao->CopyDuals){
    info = (*tao->CopyDuals)(tao,DXL,DXU,tao->data);CHKERRQ(info);
  } else {
    info=DXL->SetToZero();CHKERRQ(info);
    info=DXU->SetToZero();CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetInitialVector" 
/* ---------------------------------------------------------- */
/* @C
   TaoSetInitialVector - Sets the initial vector used by the TAO
   solver.  

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO solver
-  xx0 - vector used for initial point (set xx0==TAO_NULL for default)

   Notes:
   By default the TAO solvers use an initial starting guess of zero.  
   To provide an alternative initial guess, the user must call 
   TaoSetInitialVector() before calling TaoSolve().

   This vector will replace the variable vector set in the TaoCreate() routine.

   The user is responsible for destroying this vector.

   If this vector is TAO_NULL, TAO will use the set the previous variable 
   vector to a default starting point.

   Level: developer

.seealso: TaoSolve(), TaoGetSolution()
@ */
int TaoSetInitialVector(TAO_SOLVER tao,TaoVec *xx0)
{
  int info;
  TaoTruth flag;
  TaoVec* X;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (xx0==TAO_NULL){
    //    tao->userstart=TAO_FALSE;
  } else {
    info=TaoGetSolution(tao,&X);CHKERRQ(info);
    info=X->Compatible(xx0,&flag);CHKERRQ(info);
    if (flag == TAO_FALSE){
      SETERRQ(1,"New TAO variable vector must have identical structure as the previous variable vector.");
    }
    tao->vec_sol=xx0;
    /*    info=X->CopyFrom(xx0);CHKERRQ(info);   */
    //    tao->userstart=TAO_TRUE;
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetVariableBounds"
/*@C
   TaoGetVariableBounds - Sets the vector pointers to the vectors 
   containing the upper and lower bounds on the variables.

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Output Parameters:
+  xxll - Pointer to lower bounds on all the variables
-  xxuu - Pointer to upper bounds on all the variables

   Level: advanced

.keywords: bounds, View

@*/
int TaoGetVariableBounds(TAO_SOLVER tao,TaoVec **xxll,TaoVec **xxuu)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (xxll){
    *xxll=tao->XL;
  }
  if (xxuu){
    *xxuu=tao->XU;
  }
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetTerminationReason"
/*@
   TaoSetTerminationReason - Sets the termination reason

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  reason - one of

$    TAO_CONVERGED_ATOL (2),        (res <= atol)  
$    TAO_CONVERGED_RTOL (3),        (res/res0 <= rtol) 
$    TAO_CONVERGED_TRTOL (4),       (xdiff <= trtol) 
$    TAO_CONVERGED_MINF (5),        (f <= fmin)
$    TAO_CONVERGED_USER (6),        (user defined)

$    TAO_DIVERGED_MAXITS (-2),      (its>maxits)
$    TAO_DIVERGED_NAN (-4),         (Numerical problems)
$    TAO_DIVERGED_MAXFCN (-5),      (nfunc > maxnfuncts)
$    TAO_DIVERGED_LS_FAILURE (-6),  (line search failure)
$    TAO_DIVERGED_TR_REDUCTION (-7),
$    TAO_DIVERGED_USER (-8),        (user defined)

$    TAO_CONTINUE_ITERATING  (0)

   where
+  res - residual of optimality conditions
.  res0 - initial residual of optimality conditions
.  xdiff - current trust region size
.  f - function value
.  atol - absolute tolerance
.  rtol - relative tolerance
.  its - current iterate number
.  maxits - maximum number of iterates
.  nfunc - number of function evaluations
-  maxnfuncts - maximum number of function evaluations


   Output Parameter:

   Level: intermediate

.seealso TaoGetTerminationReason(), TaoAppSetMonitor(), TaoSetMonitor()

.keywords: convergence

@*/
int TaoSetTerminationReason(TAO_SOLVER tao,TaoTerminateReason reason)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);

  tao->reason=reason;

  TaoFunctionReturn(0);
}


/* ------------ Routines to called when destroying this application ----------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDestroyRoutine"
/*@C
   TaoSetDestroyRoutine - Sets an ADDITIONAL function that will be called when
   this object is destroyed.

   Collective on TAO_SOLVER

   Input Parameters:
+  taoapp - the TAO_SOLVER solver context
.  destroy - function pointer
-  mctx - [optional] user-defined context for private data for the 
          destroy routine (may be TAO_NULL)

   Calling sequence of destroy:
$     int mydestroy(void *mctx)

.    mctx - [optional] destroy context


   Level: intermediate

.keywords: destroy

.seealso: TaoSetMonitor(), TaoDestroy()
@*/
int TaoSetDestroyRoutine(TAO_SOLVER tao,int (*destroy)(void*),void *mctx)
{
  TaoFunctionBegin;
  if (destroy){
    if (tao->numberdestroyers >= MAX_TAO_DESTROY) {
      SETERRQ(1,"TAO ERRROR: Too many TAO destroy routines set");
    }
    
    tao->userdestroy[tao->numberdestroyers]           = destroy;
    tao->userctxdestroy[tao->numberdestroyers++]      = mctx;
  }
  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptionsPrefix"
/*@C
   TaoSetOptionsPrefix - Sets the prefix used for searching for all
   TAO options in the database.


   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
-  prefix - the prefix string to prepend to all TAO option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   For example, to distinguish between the runtime options for two
   different TAO solvers, one could call
.vb
      TaoSetOptionsPrefix(ksp1,"sys1_")
      TaoSetOptionsPrefix(ksp2,"sys2_")
.ve

   This would enable use of different options for each system, such as
.vb
      -sys1_tao_method blmvm -sys1_tao_gtol 1.e-3
      -sys2_tao_method lmvm  -sys2_tao_gtol 1.e-4
.ve


   Level: advanced

.keywords: options

.seealso: TaoAppendOptionsPrefix(), TaoGetOptionsPrefix()
@*/
int TaoSetOptionsPrefix(TAO_SOLVER tao, const char prefix[])
{
  int  info;
  TaoFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = PetscObjectSetOptionsPrefix((PetscObject)tao,prefix); CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppendOptionsPrefix"
/*@C
   TaoAppendOptionsPrefix - Appends to the prefix used for searching for all
   TAO options in the database.


   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
-  prefix - the prefix string to prepend to all TAO option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.


   Level: advanced

.keywords: options

.seealso: TaoSetOptionsPrefix(), TaoGetOptionsPrefix()
@*/
int TaoAppendOptionsPrefix(TAO_SOLVER tao, const char prefix[])
{
  int  info;
  TaoFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = PetscObjectAppendOptionsPrefix((PetscObject)tao,prefix); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoGetOptionsPrefix"
/*@C
  TaoGetOptionsPrefix - Gets the prefix used for searching for all 
  TAO options in the database

  Not Collective

  Input Parameters:
. tao - the TAO_SOLVER context
  
  Output Parameters:
. prefix - pointer to the prefix string used is returned

  Notes: On the fortran side, the user should pass in a string 'prefix' of
  sufficient length to hold the prefix.

  Level: advanced

.keywords: options

.seealso: TaoSetOptionsPrefix(), TaoAppendOptionsPrefix()
@*/
int TaoGetOptionsPrefix(TAO_SOLVER tao, const char *prefix[])
{
  int info;
  TaoFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = PetscObjectGetOptionsPrefix((PetscObject)tao,prefix); CHKERRQ(info);
  TaoFunctionReturn(0);
}
