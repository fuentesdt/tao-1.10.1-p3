#include "src/tao_impl.h"      /*I "tao_solver.h"  I*/

#undef __FUNCT__  
#define __FUNCT__ "TaoSetMethodFromOptions"
/*@
  TaoSetMethodFromOptions - Sets the TAO_SOLVER solver type from the options database, 
  or sets a default if no method has been specified.

   Collective on TAO_SOLVER

   Input Parameter:
.  solver - the TAO_SOLVER solver context

   Options Database Keys:
.  -tao_method <type> - tao_nls, tao_ntr, tao_ntl, tao_lmvm, tao_cg, tao_tron, etc.

   Level: intermediate

.keywords: method, options, database

.seealso: TaoSetOptions()
@*/
int TaoSetMethodFromOptions(TAO_SOLVER solver)
{
  char       type[256];
  int        info;
  TaoTruth flg;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  if (solver->setupcalled) SETERRQ(1,"Must call prior to TaoSetUp()");

  info = TaoOptionString("-tao_method","TaoMethod","TaoSetMethod","",type,256, &flg);CHKERRQ(info);
  if (flg) {
    info = TaoSetMethod(solver,(TaoMethod) type);CHKERRQ(info);
  }
  if (!((PetscObject)solver)->type_name) {
    info = TaoSetMethod(solver,"tao_lmvm");CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetFromOptions"
/*@
   TaoSetFromOptions - Sets many TAO_SOLVER parameters from the command line arguments.  
   This command does not set the solver type.

   Collective on TAO_SOLVER

   Input Parameter:
.  solver - the TAO_SOLVER solver context

   Options Database Keys:
+  -tao_stol - convergence tolerance in terms of the norm
                of the change in the solution between steps
.  -tao_fatol <fatol> - absolute tolerance of residual norm
.  -tao_frtol <frtol> - relative decrease in tolerance norm from initial
.  -tao_max_its <max_its> - maximum number of iterations
.  -tao_max_funcs <max_funcs> - maximum number of function evaluations
.  -tao_trtol <trtol> - trust region tolerance
.  -tao_trust0 <radius> - initial trust region radius
.  -tao_no_convergence_test - skip convergence test in minimization solver;
                               hence iterations will continue until max_it
                               or some other criterion is reached. Saves expense
                               of convergence test
.  -tao_monitor - prints residual norm at each iteration 
.  -tao_vecmonitor - plots solution at each iteration
.  -tao_vecmonitor_update - plots update to solution at each iteration 
.  -tao_xmonitor - plots residual norm at each iteration 
.  -tao_fd - use finite differences to compute Hessian; very slow, only for testing
-  -tao_mf_ksp_monitor - if using matrix-free multiply then print h at each KSP iteration

   Notes:
   To see all options, run your program with the -help option or consult
   the users manual.

   Level: developer

.keywords: options, converence, monitor, view, database

.seealso: TaoSetMethodFromOptions()
@*/
int TaoSetFromOptions(TAO_SOLVER solver)
{
  TaoTruth flg;
  int      info;
  char      type[256];

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

//  info = TaoMethodsList("-tao_method","Select TAO method","TaoSetMethod",0,type,256,0);CHKERRQ(info);
  
  info = TaoOptionString("-tao_method","TaoMethod","TaoSetMethod","",type,256, &flg);CHKERRQ(info);
  if (flg) {
    info = TaoSetMethod(solver,(TaoMethod) type);CHKERRQ(info);
  }
/*  if (!((PetscObject)solver)->type_name) {
    info = TaoSetMethod(solver,"tao_lmvm");CHKERRQ(info);
    }*/

  if (solver->setfromoptions) {
    info = (*solver->setfromoptions)(solver,solver->data); CHKERRQ(info);
  }


  info = TaoOptionName("-tao_view","view TAO_SOLVER info after each minimization has completed","TaoView",&flg);CHKERRQ(info);
  if (flg) solver->viewtao = TAO_TRUE;
  info = TaoOptionName("-tao_kspview","view the Linear Solver used by the solver after minimization has completed","TaoViewLinearSolver",&flg);CHKERRQ(info);
  if (flg) solver->viewksptao = TAO_TRUE;
  
  info = TaoOptionDouble("-tao_fatol","Stop if solution within","TaoSetTolerances",solver->fatol,&solver->fatol,&flg);CHKERRQ(info);
  info = TaoOptionDouble("-tao_frtol","Stop if relative solution within","TaoSetTolerances",solver->frtol,&solver->frtol,&flg);CHKERRQ(info);
  info = TaoOptionDouble("-tao_catol","Stop if constraints violations within","TaoSetTolerances",solver->catol,&solver->catol,&flg);CHKERRQ(info);
  info = TaoOptionDouble("-tao_crtol","Stop if relative contraint violations within","TaoSetTolerances",solver->crtol,&solver->crtol,&flg);CHKERRQ(info);
  info = TaoOptionDouble("-tao_gatol","Stop if norm of gradient less than","TaoSetGradientTolerances",solver->gatol,&solver->gatol,&flg);CHKERRQ(info);
  info = TaoOptionDouble("-tao_grtol","Stop if norm of gradient divided by the function value is less than","TaoSetGradientTolerances",solver->grtol,&solver->grtol,&flg);CHKERRQ(info); 
  info = TaoOptionDouble("-tao_gttol","Stop if the norm of the gradient is less than the norm of the initial gradient times","TaoSetGradientTolerances",solver->gttol,&solver->gttol,&flg);CHKERRQ(info); 
  info = TaoOptionInt("-tao_max_its","Stop if iteration number exceeds",
		       "TaoSetMaximumIterates",solver->max_its,&solver->max_its,
		       &flg);CHKERRQ(info);
  info = TaoOptionInt("-tao_max_funcs","Stop if number of function evaluations exceeds","TaoSetMaximumFunctionEvaluations",solver->max_funcs,&solver->max_funcs,&flg);
  info = TaoOptionDouble("-tao_fmin","Stop if function less than","TaoSetFunctionLowerBound",solver->fmin,&solver->fmin,&flg);
  info = TaoOptionDouble("-tao_steptol","Stop if step size or trust region radius less than","TaoSetTrustRegionRadius",solver->trtol,&solver->trtol,&flg);CHKERRQ(info);
  info = TaoOptionDouble("-tao_trust0","Initial trust region radius","TaoSetTrustRegionRadius",solver->trust0,&solver->trust0,&flg);CHKERRQ(info);

  /*
  info = (*PetscHelpPrintf)(solver->comm," TAO_SOLVER Monitoring Options: Choose any of the following\n");CHKERRQ(info);
  */

  info = TaoOptionName("-tao_unitstep","Always use unit step length","TaoCreateUnitLineSearch",&flg);
  if (flg){info=TaoCreateUnitLineSearch(solver);CHKERRQ(info);}
  info = TaoOptionName("-tao_lmvmh","User supplies approximate hessian for LMVM solvers","TaoLMVMSetH0",&flg);
  if (flg){info=TaoBLMVMSetH0(solver,TAO_TRUE);CHKERRQ(info);info=TaoLMVMSetH0(solver,TAO_TRUE);CHKERRQ(info);}
  
  info = TaoOptionName("-tao_view_hessian","view Hessian after each evaluation","None",&flg);CHKERRQ(info);
  if (flg) solver->viewhessian = TAO_TRUE;
  info = TaoOptionName("-tao_view_gradient","view gradient after each evaluation","None",&flg);CHKERRQ(info);
  if (flg) solver->viewgradient = TAO_TRUE;
  info = TaoOptionName("-tao_view_jacobian","view jacobian after each evaluation","None",&flg);CHKERRQ(info);  
  if (flg) solver->viewjacobian = TAO_TRUE;
  info = TaoOptionName("-tao_view_constraints","view constraint function after each evaluation","None",&flg);CHKERRQ(info);  
  if (flg) solver->viewvfunc = TAO_TRUE;

  info = TaoOptionName("-tao_cancelmonitors","cancel all monitors hardwired in code","TaoClearMonitor",&flg);CHKERRQ(info);
  if (flg) {info = TaoClearMonitor(solver);CHKERRQ(info);}
  info = TaoOptionName("-tao_monitor","Use the default convergence monitor","TaoSetMonitor",&flg);CHKERRQ(info);
  if (flg && solver->defaultmonitor) {
    info = TaoSetMonitor(solver,solver->defaultmonitor,TAO_NULL);CHKERRQ(info);
  }
  info = TaoOptionName("-tao_smonitor","Use short monitor","None",&flg);CHKERRQ(info);
  if (flg) {info = TaoSetMonitor(solver,TaoDefaultSMonitor,TAO_NULL);CHKERRQ(info);}
  info = TaoOptionName("-tao_vecmonitor","Plot solution vector at each iteration","TaoVecViewMonitor",&flg);CHKERRQ(info);
  if (flg) {info = TaoSetMonitor(solver,TaoVecViewMonitor,TAO_NULL);CHKERRQ(info);}
  info = TaoOptionName("-tao_vecmonitor_update","plots step direction at each iteration","TaoVecViewMonitorUpdate",&flg);CHKERRQ(info);
  if (flg) {info = TaoSetMonitor(solver,TaoVecViewMonitorUpdate,TAO_NULL);CHKERRQ(info);}


  //  info = PetscOptionsEnd();CHKERRQ(info);

  //  info = TaoSetLinearSolverOptions(solver);CHKERRQ(info);

  TaoFunctionReturn(0); 
}


/* -----------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoCreateFull"
/*
   TaoCreateFull - Creates a TAO_SOLVER context.

   Collective on MPI_Comm

   Input Parameters:
+  method - A TAO method.  
.  prefix - a prefix to prepend to all option names (usually TAO_NULL)
-  comm - MPI communicator

   Output Parameter:
.  newsolver - the new TAO_SOLVER context

   Options Database Keys:
.   -tao_method - select which method TAO should use

   Level: developer

.keywords: Create, solver, method, context

.seealso: TaoCreate(), TaoSolve(), TaoSetMethod(), TaoDestroy()
*/
int TaoCreateFull(TaoMethod method, const char* prefix, MPI_Comm comm, TAO_SOLVER *newsolver)
{
  TAO_SOLVER solver;
  int info;

  TaoFunctionBegin;

  *newsolver = 0;

  info = TaoInitialize(0,0,0,0);CHKERRQ(info);CHKERRQ(info);

  info = TaoObjectCreate(newsolver,comm);CHKERRQ(info);
  solver=*newsolver;

  //  info = TaoSetOptionsPrefix(solver,prefix);CHKERRQ(info);

  solver->vec_sol=0;
  solver->hessian=0;
  solver->vfunc=0;
  solver->jacobian=0;
  solver->RXL=0;
  solver->RXU=0;
  solver->CA=0;

  info = TaoResetSolver(solver); CHKERRQ(info);  /* Set some pointers to NULL */
  info = TaoSetDefaultParameters(solver);CHKERRQ(info);
  info = TaoSetDefaultStatistics(solver);CHKERRQ(info);
  info = TaoSetDefaultMonitors(solver);CHKERRQ(info);

  *newsolver = solver;
  if (method) {
    info=TaoSetMethod(solver,method);CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}


/* -----------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate"
/*@C
   TaoCreate - Creates a TAO_SOLVER context.

   Collective on MPI_Comm

   Input Parameters:
+  comm - MPI communicator
-  method - A TAO method.  

   Output Parameter:
.  newsolver - the new TAO_SOLVER context

   Options Database Keys:
.   -tao_method - select which method TAO should use

   Available methods include:
+    tao_nls - Newton's method with line search for unconstrained minimization
.    tao_ntr - Newton's method with trust region for unconstrained minimization
.    tao_ntl - Newton's method with trust region, line search for unconstrained minimization
.    tao_lmvm - Limited memory variable metric method for unconstrained minimization
.    tao_cg - Nonlinear conjugate gradient method for unconstrained minimization
.    tao_nm - Nelder-Mead algorithm for derivate-free unconstrained minimization
.    tao_tron - Newton Trust Region method for bound constrained minimization
.    tao_gpcg - Newton Trust Region method for quadratic bound constrained minimization
.    tao_blmvm - Limited memory variable metric method for bound constrained minimization
.    tao_kt - Formulate a bound constrained problem as a complementarity problem
.    tao_bqpip - Interior point method for quadratic bound constrained minimization
.    tao_ssils - Infeasible semismooth method with a linesearch for complementarity problems
-    tao_ssfls - Feasible semismooth method with a linesearch for complementarity problems

   Level: beginner

   Note: 
   If the second argument specifies a TaoMethod, quotation marks should 
   surround the method.

   Note:
   The TaoMethod can be TAO_NULL (C/C++) or 
   TAO_NULL_CHARACTER (Fortran), in which case the
   method will be specified by the runtime option -tao_method

   If a particular optimization method is specified at runtime by
   the option '-tao_method', this choice will be used instead of
   any default that may have been specified as the input parameter
   "method" to this routine.

.keywords: Create, solver, method, context

.seealso: TaoSolve(), TaoSetMethod(), TaoSetApplication(), TaoDestroy()
@*/
int TaoCreate(MPI_Comm comm, TaoMethod method, TAO_SOLVER *newsolver)
{
  int info;

  TaoFunctionBegin;

  *newsolver = 0;

  info = TaoCreateFull(method,0,comm,newsolver);CHKERRQ(info);
  info = TaoSetMethodFromOptions(*newsolver);CHKERRQ(info);

  TaoFunctionReturn(0);
}



/* ----- Routines to initialize and destroy a minimization solver ---- */

#undef __FUNCT__  
#define __FUNCT__ "TaoSetDefaultParameters"
/*@
   TaoSetDefaultParameters - Set the parameters used by all TAO solvers to a default value.  These
   parameter include convergence tolerances.  This routine is called before setting the
   method used by TAO

   Collective on TAO_SOLVER

   Input Parameters:
.  solver - the TAO_SOLVER solver context

   Level: developer

.keywords: options, defaults

.seealso: TaoCreate(), TaoSetDefaultStatistics(), TaoSetDefaultMonitors()
@*/
int TaoSetDefaultParameters(TAO_SOLVER solver){

  TaoFunctionBegin;
  solver->max_its            = 0;
  solver->max_funcs	     = 100000000;
  solver->fatol              = 0.0;
  solver->frtol              = 0.0;
  solver->catol              = 0.0;
  solver->crtol              = 0.0;
  solver->gatol              = 0.0;
  solver->grtol              = 0.0;
  solver->xtol		     = 0.0;
  solver->trtol              = 0.0;
  solver->fmin               = -1.e30;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetDefaultStatistics"
/*@
   TaoSetDefaultStatistics - Initialize the statistics used by TAO for all of the solvers.
   These statistics include the iteration number, residual norms, and convergence status.
   This routine gets called before solving each optimization problem.

   Collective on TAO_SOLVER

   Input Parameters:
.  solver - the TAO_SOLVER solver context

   Level: developer

.keywords: options, defaults

.seealso: TaoCreate(), TaoSetDefaultParameters(), TaoSetDefaultMonitors(), TaoSolve()
@*/
int TaoSetDefaultStatistics(TAO_SOLVER solver){

  TaoFunctionBegin;
  solver->iter               = 0;
  solver->fc                 = 0;
  solver->norm		     = 0.0;
  solver->norm0		     = 0.0;
  solver->cnorm		     = 0.0;
  solver->cnorm0	     = 0.0;
  solver->nfuncs             = 0;
  solver->ngrads             = 0;
  solver->nfgrads            = 0;
  solver->nhesss             = 0;
  solver->nvfunc             = 0;
  solver->njac               = 0;
  solver->linear_its         = 0;
  solver->lsflag             = 0;
  solver->reason             = TAO_CONTINUE_ITERATING;
  solver->step               = 1.0e+30;
  if (solver->conv_hist_reset == TAO_TRUE) solver->conv_hist_len = 0;
  TaoFunctionReturn(0);

}
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDefaultMonitors"
/*@
   TaoSetDefaultMonitors - Set the default monitors and viewing options available in TAO.
   This routine is generally called only in TaoCreate().

   Collective on TAO_SOLVER

   Input Parameters:
.  solver - the TAO_SOLVER solver context

   Level: developer

.keywords: options, defaults

.seealso: TaoCreate(), TaoSetDefaultStatistics(), TaoSetDefaultParameters(), TaoSetMonitor(), TaoView()
@*/
int TaoSetDefaultMonitors(TAO_SOLVER solver){

  TaoFunctionBegin;
  solver->numbermonitors     = 0;
  solver->viewhessian        = TAO_FALSE;
  solver->viewgradient       = TAO_FALSE;
  solver->viewjacobian       = TAO_FALSE;
  solver->viewvfunc          = TAO_FALSE;
  solver->viewtao            = TAO_FALSE;
  solver->viewksptao        = TAO_FALSE;

  solver->converged          = TaoConverged_Default;
  solver->defaultmonitor     = TaoDefaultMonitor;

  solver->conv_hist_len      = 0;
  solver->conv_hist_max      = 0;
  solver->conv_hist          = TAO_NULL;
  solver->conv_hist_its      = TAO_NULL;
  solver->conv_hist_reset    = TAO_TRUE;

  solver->numberdestroyers   =0;

  TaoFunctionReturn(0);

}

/* ----- Routines to initialize and destroy a minimization solver ---- */

#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp"
/*@
   TaoSetUp - Sets up the internal data structures for the later use
   of a minimization solver.

   Collective on TAO_SOLVER

   Input Parameters:
.  solver - the TAO_SOLVER solver context

   Notes:
   For basic use of the TAO_SOLVER solvers the user need not explicitly call
   TaoSetUp(), since these actions will automatically occur during
   the call to TaoSolve().  However, if one wishes to control this
   phase separately, TaoSetUp() should be called after TaoCreate()
   and optional routines of the form TaoSetXXX(), but before TaoSolve().  

   Level: developer

.keywords: Solve, setup

.seealso: TaoCreate(), TaoSolve(), TaoSetDown(), TaoDestroy()
@*/
int TaoSetUp(TAO_SOLVER solver)
{
  int      info;
  TaoTruth flag;
  TaoVec   *xx,*dx;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  //  info = TaoSetOptions(solver);CHKERRQ(info);
  if (!solver->set_method_called) {
    SETERRQ(1,"Must explicitly call TaoSetMethod() or TaoSetMethodFromOptions() before TaoSolve()");
  }

  /* Determine if the solver has already been set up with structures of the right dimension */
  if ( solver->setupcalled==TAO_TRUE) {
    info = TaoGetSolution(solver,&xx);CHKERRQ(info);
    if (!xx){SETERRQ(1,"Must explicitly call TaoSetApplication() and Set Variable Vector");}
    info = TaoGetStepDirectionVector(solver,&dx);CHKERRQ(info);
    if (dx){
      info = xx->Compatible(dx,&flag); CHKERRQ(info);
    } else {
      flag=TAO_FALSE;
    }
    if (flag==TAO_TRUE){
      info = TaoGetGradient(solver,&dx);CHKERRQ(info);
      if (dx){
	info = xx->Compatible(dx,&flag); CHKERRQ(info);
      } 
    }
    if (flag==TAO_FALSE){  /* Setup done, but data structures of wrong size */
      info = TaoSetDown(solver); CHKERRQ(info);
    }
  }

  if ( solver->setupcalled==TAO_FALSE) {
    if (solver->setup) {
      info = (*solver->setup)(solver,solver->data);CHKERRQ(info);    
    }
  } 
  solver->setupcalled=TAO_TRUE;
  info = TaoLineSearchSetUp(solver);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy"
/*@
   TaoDestroy - Destroys the TAO solver that was created with TaoCreate().

   Collective on TAO_SOLVER

   Input Parameter:
.  solver - the TAO_SOLVER solver context

   Level: beginner

.keywords: Destroy

.seealso: TaoCreate(), TaoSolve()
@*/
int TaoDestroy(TAO_SOLVER solver)
{
  int i,info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  if (--((PetscObject)solver)->refct > 0) TaoFunctionReturn(0);

  for (i=0; i< solver->numberdestroyers; i++){
    info = (*solver->userdestroy[i])(solver->userctxdestroy[i]); CHKERRQ(info);
  }

  info = TaoResetSolver(solver);CHKERRQ(info);

  info = TaoObjectDestroy(solver); CHKERRQ(info);
  
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown"
/*@
   TaoSetDown - Take down the data structures created in TaoSetUp().
   These structures typically include the work vectors, and linear solver.

   Collective on TAO_SOLVER

   Input Parameter:
.  solver - the TAO_SOLVER solver context

   Level: advanced

.keywords: Destroy

.seealso: TaoSetUp(), TaoDestroy()
@*/
int TaoSetDown(TAO_SOLVER solver)
{
  int info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  if (solver->setupcalled){
    if (solver->setdown) {info = (*(solver)->setdown)(solver,solver->data);CHKERRQ(info);}
  }
  info = TaoSetLagrangianGradientVector(solver,0);CHKERRQ(info);
  info = TaoSetStepDirectionVector(solver,0);CHKERRQ(info);
  info = TaoSetVariableBounds(solver,0,0);CHKERRQ(info);
  solver->ksp = TAO_NULL;
  solver->setupcalled=TAO_FALSE;

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoResetSolver"
/*@
   TaoResetSolver - Take down the data structures created in TaoCreate__XXX().
   This routine destroys the line search, and the solver context.  It
   also set many of the solver routines for solving, options, duals, viewing,
   setup, and destroy to TAO_NULL.

   Collective on TAO_SOLVER

   Input Parameter:
.  solver - the TAO_SOLVER solver context

   Level: advanced

.keywords: Destroy

.seealso: TaoCreate(), TaoSetMethod(), TaoSetDown(), TaoDestroy()
@*/
int TaoResetSolver(TAO_SOLVER solver)
{
  int info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  info = TaoSetDown(solver); CHKERRQ(info);
  if (solver->data){
    info = TaoFree(solver->data);CHKERRQ(info);
  }
  solver->data = 0;
  solver->CopyDuals=0;
  solver->solve=0;
  solver->data=0;
  solver->view=0;
  solver->setup=0;
  solver->setdown=0;
  solver->setfromoptions=0;

  info = TaoLineSearchDestroy(solver);CHKERRQ(info);
  info = TaoMeritFunctionDestroy(solver);CHKERRQ(info);

  /* Set Default Parameters */
  info = TaoSetDefaultParameters(solver); CHKERRQ(info);
  info = TaoSetDefaultMeritFunction(solver); CHKERRQ(info);

  solver->set_method_called  = TAO_FALSE;

  TaoFunctionReturn(0);
}



/* --------- Internal routines for TAO_SOLVER Package --------- */

#undef __FUNCT__  
#define __FUNCT__ "TaoSetMethod"
/*@C
   TaoSetMethod - Sets the method for the unconstrained minimization solver.  

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  type - a known method

   Options Database Key:
.  -tao_method <type> - Sets the method; use -help for a list
   of available methods (for instance, "-tao_method tao_lmvm" or 
   "-tao_method tao_tron")

   Available methods include:
+    tao_nls - Newton's method with line search for unconstrained minimization
.    tao_ntr - Newton's method with trust region for unconstrained minimization
.    tao_ntl - Newton's method with trust region, line search for unconstrained minimization
.    tao_lmvm - Limited memory variable metric method for unconstrained minimization
.    tao_cg - Nonlinear conjugate gradient method for unconstrained minimization
.    tao_nm - Nelder-Mead algorithm for derivate-free unconstrained minimization
.    tao_tron - Newton Trust Region method for bound constrained minimization
.    tao_gpcg - Newton Trust Region method for quadratic bound constrained minimization
.    tao_blmvm - Limited memory variable metric method for bound constrained minimization
.    tao_kt - Formulate a bound constrained problem as a complementarity problem
.    tao_bqpip - Interior point method for quadratic bound constrained minimization
.    tao_ssils - Infeasible semismooth method with a linesearch for complementarity problems
-    tao_ssfls - Feasible semismooth method with a linesearch for complementarity problems

  Level: intermediate

.keywords: method, Create, solve

.seealso: TaoCreate(), TaoGetMethod()

@*/
int TaoSetMethod(TAO_SOLVER solver,TaoMethod type)
{
  int info;
  int (*r)(TAO_SOLVER);
  TaoTruth issame;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  info = TaoCompareMethod(solver,type,&issame);CHKERRQ(info);
  if (issame) TaoFunctionReturn(0);
  info = TaoResetSolver(solver); CHKERRQ(info);
  info = TaoFindSolver(solver,type,&r);CHKERRQ(info);  

  if (!r) SETERRQ1(1,"Unable to find requested TAO_SOLVER type %s",type);
  info = (*r)(solver);CHKERRQ(info);

  solver->set_method_called = TAO_TRUE;

  TaoFunctionReturn(0); 
}


/* --------------------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetTaoDualVariablesRoutine"
/*@C
  TaoSetTaoDualVariablesRoutine - Set a routine that can be called
  to compute the dual variables on the lower and upper bounds of the
  variables.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  dual - Dual variables routine

   Note: The calling sequence of the dual routine passes
   the TAO_SOLVER object in the first argument, an (optional)
   pointer to a TaoVec to put the duals of the lower bounds,
   an (optional) pointer to a TaoVec to put the duals of 
   the upper bounds, and the solver context passed to TAO
   in TaoSetSolver().

   Level: developer

.keywords: TAO_SOLVER, duals

.seealso: TaoSetTaoSolveRoutine()
@*/
int TaoSetTaoDualVariablesRoutine(TAO_SOLVER tao, 
				  int (*duals)(TAO_SOLVER,TaoVec*,TaoVec*,void*))
{
  TaoFunctionBegin;
  tao->CopyDuals  = duals;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetTaoSolveRoutine"
/*@C
   TaoSetTaoSolveRoutine - Sets the routine that will solve an optimization application

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER context
.  solve - routine that applies the algorithm
-  data - solver data structure (optional)

   Level: developer

   Note:
   This routine is generally used within a "TaoCreate_XXX" routine.
   TAO will call this routine as part of the TaoSolve() routine.

   Note:
   The first and third arguments of this routine will be used as
   the arguments of the solver routine provided here.

.keywords: TAO_SOLVER, solve

.seealso: TaoCreate(), TaoSolve(), TaoGetSolverContext()
@*/
int TaoSetTaoSolveRoutine(TAO_SOLVER tao, 
			  int (*solve)(TAO_SOLVER,void*), void*ctx)
{  
  TaoFunctionBegin;
  tao->solve		  = solve;
  tao->data	          = ctx;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetTaoSetUpDownRoutines"
/*@C
   TaoSetTaoSetUpDownRoutines - Sets the routines that setup and destroy solver data structures

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
.  setup - routine that creates the work vectors in a solver.
-  setdown - the routine that will destroy the work vectors of a solver

   Note:
   This routine is generally called within a "TaoCreate_XXX" routine.
   The routines set here will be called in the TaoSetApplication() and
   TaoDestroy() routines, respectively.  Vectors and other data structures
   needed by the solver can be created and destroyed within the TaoSolve_XXX()
   routine, or before and after this routine.  The advantage to doing it
   before and after is that the solver can be called multiple times
   without reallocated these structures -- improving efficiency.

   Note:
   When the 'setup' routine is called, the solution vector, and other
   data will be available to clone.
   
   Note:
   When TAO calls these routines, the second arguement will be the
   context specified in TaoSetTaoSolveRoutine().

   Level: developer

.keywords: TAO_SOLVER, setup, destroy

.seealso: TaoCreate(), TaoSetUp(), TaoSetDown(), TaoDestroy(), TaoSetTaoSolveRoutine()
@*/
int TaoSetTaoSetUpDownRoutines(TAO_SOLVER tao, 
				 int (*setup)(TAO_SOLVER,void*),
				 int (*setdown)(TAO_SOLVER,void*))
{
  TaoFunctionBegin;
  tao->setup		  = setup;
  tao->setdown		  = setdown;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetTaoViewRoutine"
/*@C
   TaoSetTaoViewRoutine - Sets the routine that will display information
   about the optimization solver

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  view - routine that views the solver

   Note:
   This routine is generally used within a "TaoCreate_XXX" routine.
   TAO will call this routine as part of the TaoView() routine.

   Note:
   When TAO calls these routines, the second arguement will be the
   context specified in TaoSetTaoSolveRoutine().

   Level: developer

.keywords: TAO_SOLVER, view

.seealso: TaoCreate(), TaoView(), TaoSetTaoSolveRoutine()
@*/
int TaoSetTaoViewRoutine(TAO_SOLVER tao, 
			 int (*view)(TAO_SOLVER,void*))
{  
  TaoFunctionBegin;
  tao->view               = view;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetTaoOptionsRoutine"
/*@C
   TaoSetTaoSolveRoutine - Sets the routine that will set the options for
   the optimization solver.

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER context
-  options - routine that views the solver

   Note:
   This routine is generally used within a "TaoCreate_XXX" routine.
   TAO will call this routine as part of the TaoSetOptions() routine.

   Note:
   When TAO calls these routines, the second argument will be the
   context specified in TaoSetTaoSolveRoutine().

   Level: developer

.keywords: TAO_SOLVER, solve

.seealso: TaoCreate(), TaoSetOptions(), TaoSetTaoSolveRoutine().
@*/
int TaoSetTaoOptionsRoutine(TAO_SOLVER tao, 
			    int (*options)(TAO_SOLVER,void*))
{  
  TaoFunctionBegin;
  tao->setfromoptions     = options;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetApplication"
/*@C
   TaoGetApplication - Gets the user defined context for 
   the minimization solvers.  

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
-  taoapp - user application context

   Level: advanced

.keywords: application, context

@*/
int TaoGetApplication(TAO_SOLVER tao, TaoApplication **taoapp){
  TaoFunctionBegin;
  *taoapp=tao->taoappl;
  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoGetSolverContext"
/*@C
   TaoGetSolverContext - Gets the solver context for 
   the minimization solvers.  

   Collective on TAO_SOLVER

   Input Parameter:
+  tao - the TAO_SOLVER solver context
-  type - the name of the method

   Output Parameter:
.  solverctx - solver context IF the type matches.

   Level: developer

.keywords: application, context

.seealso: TaoSetTaoSolveRoutine()

@*/
int TaoGetSolverContext(TAO_SOLVER tao, TaoMethod type, void **solverctx){
  int info;
  TaoTruth issame;
  TaoFunctionBegin;
  if (solverctx){
    info = TaoCompareMethod(tao,type,&issame);CHKERRQ(info);
    if (issame) *solverctx=tao->data;
    else        *solverctx=0;
  }
  TaoFunctionReturn(0);
}

