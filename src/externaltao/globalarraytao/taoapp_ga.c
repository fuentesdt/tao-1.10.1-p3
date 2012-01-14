//File: taoapp_ga.c

/**************************************************************

Author: Limin Zhang, Ph.D.
        Mathematics Department
        Columbia Basin College
        Pasco, WA 99301
        Limin.Zhang@cbc2.org

Mentor: Jarek Naplocha, Ph.D.
        Environmental Molecular Science Laboratory
        Pacific Northwest National Laboratory
        Richland, WA 99352

Date: 7/22/2002

Purpose:
      to design and implement an application interface between TAO and global arrays.
**************************************************************/

#include "taoapp_ga.h"         /*I "taoapp_ga.h" I*/




#undef __FUNCT__ 
#define __FUNCT__ "TaoGAApplication::TaoGAApplication"
TaoGAApplication::TaoGAApplication(MPI_Comm mpicomm){

  this->comm=mpicomm;
  this->taox=0;
  this->V = 0;
  this->usrfctx=0; 
  this->usrgctx=0; 
  this->usrfgctx=0; 
  this->HesMat=0; 
  this->taoh=0;
  this->usrhctx=0;
  this->numbermonitors = 0;

  this->computeumfunction=0; 
  this->computegradient=0;
  this->computefunctiongradient=0; 
  this->computehessian=0;
  this->computebounds=0;
  for (int i=0;i<MAX_TAO_MONITORS;i++) {
    this->monitor[i] = 0;
    this->monitorcontext[i] = 0;
  }

  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::~TaoGAApplication"
TaoGAApplication::~TaoGAApplication(){
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::GetVariableVector"
int TaoGAApplication::GetVariableVector(TaoVec **xx){
  int info;
  TaoFunctionBegin;
  if (this->V == 0) {
    SETERRQ(1,"Variable vector is not set. Did you call TaoGAAppSetInitialSolutionVec?");
  }
  if (this->taox == 0) {
    info = TaoWrapGaVec(this->V,&this->taox);
  }
  *xx= this->taox;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::GetHessianMatrix"
int TaoGAApplication::GetHessianMatrix(TaoMat **HH){
  TaoFunctionBegin;
  *HH = this->taoh;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::EvaluateVariableBounds"
int TaoGAApplication::EvaluateVariableBounds(TaoVec *xxll, TaoVec *xxuu) {
  int info;
  GAVec XL=0, XU=0;
  TaoFunctionBegin;
  if (this->computebounds) {
    info = TaoVecGetGaVec(xxll, &XL); CHKERRQ(info);
    info = TaoVecGetGaVec(xxuu, &XU); CHKERRQ(info);
    info = (*this->computebounds)(this, XL, XU, this->usrvbctx); CHKERRQ(info);
  } 
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::EvaluateObjectiveFunction"
int TaoGAApplication::EvaluateObjectiveFunction(TaoVec *xx, double *ff){
  int     info;
  GAVec X = 0;
  TaoVecGetGaVec(xx, &X);
  GAVec G;
  TaoVec *gg;
  
  TaoFunctionBegin;

  if (this->computeumfunction){
    info = (*this->computeumfunction)(this,X,ff,this->usrfctx); CHKERRQ(info);
  } else if (this->computefunctiongradient) {
    info = xx->Clone(&gg); CHKERRQ(info);
    info = TaoVecGetGaVec(gg, &G); CHKERRQ(info);
    info = (*this->computefunctiongradient)(this,X,ff,G,this->usrfgctx); CHKERRQ(info);
    info=TaoVecDestroy(gg); CHKERRQ(info);
  }
  else {
    SETERRQ(1,"function evaluation routine is not set");
  }

  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::EvaluateGradient"
int TaoGAApplication::EvaluateGradient(TaoVec *xx, TaoVec *gg){
  int     info;
  GAVec X = 0, G = 0;
  TaoVecGetGaVec(xx, &X);
  TaoVecGetGaVec(gg, &G);
  double ff;

  TaoFunctionBegin;

  if (this->computegradient){
    info = (*this->computegradient)(this,X,G,this->usrgctx); CHKERRQ(info);
  } else if ( this->computefunctiongradient ) {
    info = (*this->computefunctiongradient)(this,X,&ff,G,this->usrfgctx);
    CHKERRQ(info);
  }
  else {
    SETERRQ(1,"Gradient evaluation routine is not set.");
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::EvaluateObjectiveAndGradient"
int TaoGAApplication::EvaluateObjectiveAndGradient(TaoVec *xx, double *ff, TaoVec *gg){
  int     info;
  GAVec X = 0, G = 0;
  TaoVecGetGaVec(xx, &X);
  TaoVecGetGaVec(gg, &G);
  TaoFunctionBegin;

  if (this->computefunctiongradient){
    info = (*this->computefunctiongradient)(this,X,ff,G,this->usrfgctx);
    CHKERRQ(info);
  } else if ( this->computeumfunction && this->computegradient ) {
    info = (*this->computeumfunction)(this,X,ff,this->usrfctx); CHKERRQ(info);
    info = (*this->computegradient)(this,X,G,this->usrgctx); CHKERRQ(info);
  } else {
    SETERRQ(1,"Function and Gradient evaluation routines not set.");
  }
  

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::EvaluateHessian"
int TaoGAApplication::EvaluateHessian(TaoVec *xx, TaoMat *HH){
  int     info;

  GAVec X = 0;
  TaoVecGetGaVec(xx, &X);
  TaoMatGa *MM=(TaoMatGa*)HH;
  GAMat H=MM->GetMat(), Hpre=MM->pm_pre;;

  TaoFunctionBegin;
  if (this->computehessian){
    info = (*this->computehessian)(this,X,H,this->usrhctx); CHKERRQ(info);
    MM->pm=H;
    MM->pm_pre=Hpre;
  } else {
    SETERRQ(1,"Hessian calculation routine was not set");
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoGAApplication::InitializeVariables"
int TaoGAApplication::InitializeVariables(TaoVec *xx){
  GAVec XX0=0;

  TaoFunctionBegin; 
  if (xx){
    TaoVecGetGaVec(xx, &XX0);
  }
  if (XX0 && this->V && XX0!=this->V){
     GA_Copy(this->V,XX0); 
  }
  TaoFunctionReturn(0);
}




// C stubs

#undef __FUNCT__
#define __FUNCT__ "TaoGAApplicationCreate"
/*@C
  TaoGAApplicationCreate - Creates a TaoApplication that
uses GA data structures.   The vectors used for gradient
and other information can be a GA Vec.  The routines
for function evaluation, and derivative information can
also used GA arguments.

   Input Parameters:
.  comm - an MPI communiicator

   Output Parameters:
.  newapp - the TaoApplication structure

.seealso TaoGAAppSetObjectiveAndGradientRoutine(), TaoGAAppDestroy()

   Level: beginner

.keywords: GAApplication
@*/
int TaoGAApplicationCreate(MPI_Comm comm, TAO_GA_APPLICATION* newapp){
  *newapp=new TaoGAApplication(comm);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoGAAppDestroy"
/*@C
  TaoGAAppDestroy - Destroy the GA application 
and all of the vectors and matrices associated wit it.

   Input Parameters:
.  gaapp - the GAApplication structure

.seealso TaoGAApplicationCreate()

   Level: beginner

.keywords: Concepts: GAApplication, Destroy
@*/
int TaoGAAppDestroy(TAO_GA_APPLICATION gaapp) {
  TaoFunctionBegin;
  delete gaapp;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoGAAppSetGradientRoutine"
/*@C
   TaoGAAppSetGradientRoutine - Sets the gradient evaluation routine and gradient
   vector for use by the TAO_GA_APPLICATION routines.

   Collective on TAO_GA_APPLICATION

   Input Parameters:
+  gaapp - the TAO_GA_APPLICATION context
.  grad - gradient evaluation routine
-  ctx - [optional] user-defined function context 

   Calling sequence of func:
$     grad (TAO_GA_APPLICATION gaapp,GAVec x,GAVec g,void *ctx);

+  gaapp - the TAO_GA_APPLICATION application context
.  x - input vector
.  g - gradient vector
-  ctx - user-defined function gradient context

   Note:
   This routine should be called before TaoSetApplication()

   Level: beginner

   Options Database Keys:
.  -tao_view_gradient - view the gradient after each evaluation

.keywords: GAApplication, set, function

.seealso: TaoGAAppSetObjectiveAndGradientRoutine(), TaoGAAppSetHessianRoutine()

@*/
int TaoGAAppSetGradientRoutine(TAO_GA_APPLICATION gaapp, int (*grad)(TAO_GA_APPLICATION,GAVec,GAVec,void*),void *ctx){

  TaoFunctionBegin;
  gaapp->computegradient=grad;
  gaapp->usrgctx=ctx;
  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "TaoSetupGAApplicationSolver"
/*@
   TaoSetupGAApplicationSolver - This routine creates the vectors,
   matrices, linear solvers, and other data structures used in
   the during the optimization process.  The application provides
   the solver with an objective function, constraints, derivative 
   information, and application data structures.  These structures
   include a vector of variables, and Hessian matrix.

   Collective on TAO_SOLVER

   Input Parameters:
+  myapp - user application context
-  tao - the TAO_SOLVER solver context

   Level: intermediate

   Note: This routine should be called after TaoGAAppSetInitialSolutionVec() 

   Note: 
   This method is called during  TaoSetOptions() and TaoApplicationSolve();
   
.keywords: GAApplication

.seealso: TaoSolveGAApplication()

@*/
int TaoSetupGAApplicationSolver(TAO_GA_APPLICATION myapp, TAO_SOLVER tao ){
  int info;
  TaoFunctionBegin;
  info = TaoSetApplication(tao,myapp);CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__ 
#define __FUNCT__ "TaoSolveGAApplication"
/*@
  TaoSolveGAApplication - Find a solution to the application using a TAO solver.

   Collective on TAO_GA_APPLICATION

   Input Parameters:
.  tao - the TAO_GA_APPLICATION context

   Level: beginner

.keywords: solve

.seealso: TaoSolve();

@*/
int TaoSolveGAApplication(TAO_GA_APPLICATION gaapp, TAO_SOLVER tao){
  int info;

  TaoFunctionBegin;
  info = TaoSetupGAApplicationSolver(gaapp, tao); CHKERRQ(info);
  info = TaoSolve(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "TaoGAAppSetObjectiveAndGradientRoutine"
/*@C
   TaoGAAppSetObjectiveAndGradientRoutine - Sets a routine for function and gradient evaluation.

   Collective on TAO_GA_APPLICATION

   Input Parameters:
+  gaapp - the TAO_GA_APPLICATION context
.  funcgrad - routine for evaluating the function and gradient
-  ctx - optional user-defined context for private data for the 
         function and gradient evaluation routine (may be TAO_NULL)

   Calling sequence of funcgrad:
$     funcgrad (TAO_GA_APPLICATION tao,GAVec x,double *f,GAVec g,void *ctx);

+  tao - TAO_GA_APPLICATION application context
.  x - input vector
.  f - function value
.  g - gradient vector
-  ctx - optional user-defined context 

   Notes:
   The user may call TaoGAAppSetObjectiveAndGradientRoutine() to set a routine
   that evaluates both the function and gradient.  Alternatively, the
   user may call both TaoGAAppSetObjectiveRoutine() and TaoGAAppSetGradientRoutine() to set
   separate routines for function and gradient evaluation.  

   Using a single routine to compute the function and gradient, as
   specified via TaoGAAppSetObjectiveAndGradientRoutine(), may enable better performance
   for applications in which many of the function and gradient computations
   are identical.

   Level: beginner

   Options Database Keys:
.   -tao_view_gradient - view the gradient after each iteration

.keywords: GAApplication, set, function

.seealso: TaoGAAppSetGradientRoutine(), TaoGAAppSetObjectiveRoutine(), TaoComputeFunctionGradient()

@*/
int TaoGAAppSetObjectiveAndGradientRoutine(TAO_GA_APPLICATION gaapp, int (*funcgrad)(TAO_GA_APPLICATION,GAVec,double*,GAVec, void*),void *ctx){
  TaoFunctionBegin;
  gaapp->computefunctiongradient=funcgrad;
  gaapp->usrfgctx=ctx;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoGAAppSetObjectiveRoutine"
/*@C
   TaoGAAppSetObjectiveRoutine - Sets a routine for function evaluations.

   Collective on TAO_GA_APPLICATION

   Input Parameters:
+  gaapp - the TAO_GA_APPLICATION context
.  func - routine for evaluating the function
-  ctx - optional user-defined context for private data for the 
         function evaluation routine (may be TAO_NULL)

   Calling sequence of funcgrad:
$     func (TAO_GA_APPLICATION tao,GAVec x,double *f,void *ctx);

+  tao - TAO_GA_APPLICATION application context
.  x - input vector
.  f - function value
-  ctx - optional user-defined context 

   Notes:
   The user may call TaoGAAppSetObjectiveAndGradientRoutine() to set a routine
   that evaluates both the function and gradient.  Alternatively, the
   user may call both TaoGAAppSetObjectiveRoutine() and TaoGAAppSetGradientRoutine() to set
   separate routines for function and gradient evaluation.  

   Using a single routine to compute the function and gradient, as
   specified via TaoGAAppSetObjectiveAndGradientRoutine(), may enable better performance
   for applications in which many of the function and gradient computations
   are identical.

   Level: beginner

   Options Database Keys:
.   -tao_view_gradient - view the gradient after each iteration

.keywords: TAO_GA_APPLICATION, set, function

.seealso: TaoGAAppSetGradientRoutine(), TaoGAAppSetObjectiveAndGradientRoutine(), TaoComputeFunctionGradient()

@*/
int TaoGAAppSetObjectiveRoutine(TAO_GA_APPLICATION gaapp, int (*func)(TAO_GA_APPLICATION,GAVec,double*,void*), void *ctx) {
  TaoFunctionBegin;
  gaapp->computeumfunction = func;
  gaapp->usrfctx=ctx;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoGAAppSetHessianRoutine"
/*@C
   TaoGAAppSetHessianRoutine - Sets the function to compute the Hessian as well as the
   location to store the matrix.

   Collective on TAO_GA_APPLICATION and Mat

   Input Parameters:
+  gaapp - the TAO_GA_APPLICATION context
.  hess - Hessian evaluation routine
-  ctx - [optional] user-defined context for private data for the 
         Hessian evaluation routine (may be TAO_NULL)

   Calling sequence of hess:
$    hess (TAO_GA_APPLICATION gaapp,GAVec x,GAMat H,void *ctx);

+  gaapp - the TAO_GA_APPLICATION application context
.  x - input vector
.  H - Hessian matrix
-  ctx - [optional] user-defined Hessian context

   Options Database Keys:
.  -tao_view_hessian - view the hessian after each evaluation

   Level: beginner

.keywords: GAApplication, Hessian

.seealso: TaoGAAppSetObjectiveRoutine(), TaoGAAppSetGradientRoutine()
@*/
int TaoGAAppSetHessianRoutine(TAO_GA_APPLICATION gaapp, int (*hess)(TAO_GA_APPLICATION, GAVec, GAMat,void*),void *ctx){
  TaoFunctionBegin;
  gaapp->computehessian=hess;
  gaapp->usrhctx=ctx;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoGAAppSetVariableBoundsRoutine"
/*@C
   TaoGAAppSetVariableBoundsRoutine - Set bounds on the variables.

   Collective on TAO_GA_APPLICATION

   Input Parameters:
+  gaapp - the TAO_GA_APPLICATION context
.  func - the variable bounds evaluation routine
-  ctx - user-defined context for private data for the variable bounds
         evaluation routine (may be TAO_NULL)

   calling sequence of func:
$    func (TAO_GA_APPLICATION gaapp, GAVec xl, GAVec xu, void *ctx);

+  gaapp - the TAO_GA_APPLICATION application context
.  xl - lower bounds vector
.  xu - upper bounds vector
-  ctx - user-defined variable bounds context

.keywords: GAApplication, bounds

   Level: beginner

@*/
int TaoGAAppSetVariableBoundsRoutine(TAO_GA_APPLICATION gaapp, 
				       int (*func)(TAO_GA_APPLICATION, GAVec, GAVec, void *), void *ctx) {
  TaoFunctionBegin;
  gaapp->computebounds = func;
  gaapp->usrvbctx = ctx;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoGAAppSetInitialSolutionVec"
/*@
   TaoGAAppSetInitialSolutionVec - Sets the vector representing the variables
   and an initial guess.

   Collective on TAO_GA_APPLICATION

   Input Parameters:
+  taoapp - the TAO_GA_APPLICATION context
.  xx - variable vector that stores the solution

   Level: beginner

   Note: This vector will be used by the solver, so do not use it
   for other purposes.  The user should destroy this vector after
   solving the application.

   Note:  If the user is unaware of a decent initial solution,
   the vector should be set to zero.

   Note:  The TAO solvers will not use the contents of this 
   GAVec until the TaoSolveGAApplication() is called.  Therefore the user
   may compute an initial solution in this vector after this
   routine -- but before TaoSolveGAApplication().

.seealso:  TaoAppGetSolutionVec(), TaoAppSetObjectiveRoutine()
@*/
int TaoGAAppSetInitialSolutionVec(TAO_GA_APPLICATION gaapp, GAVec xx){
  TaoFunctionBegin;
  if (gaapp->V){
     GA_Destroy(gaapp->V);
  }  
  gaapp->V=xx;
  TaoFunctionReturn(0);
}


/* -------------- Routines to set performance monitoring options --------------- */

#undef __FUNCT__
#define __FUNCT__ "TaoGAAppSetMonitor"
/*@C
   TaoGAAppSetMonitor - Sets an ADDITIONAL function that is to be used at every
   iteration of the solver to display the interation's progress.

   Collective on TAO_GA_APPLICATION

   Input Parameters:
+  gaapp - the TAO_GA_APPLICATION application context
.  mymonitor - monitoring routine
-  mctx - [optional] user-defined context for private data for the
          monitor routine (may be TAO_NULL)

   Calling sequence of mymonitor:
$     int mymonitor(TAO_DA_APPLICATION gaapp, void *mctx)

+    gaapp - the TAO_GA_APPLICATION application context
-    mctx - [optional] monitoring context

   Notes:
   Several different monitoring routines may be set by calling
   TaoGAAppSetMonitor() mutiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.keywords: GAApplication, monitor, View

.seealso: TaoGAAppMonitor()
@*/
int TaoGAAppSetMonitor(TAO_GA_APPLICATION gaapp, int (*mymonitor)(TAO_GA_APPLICATION, void*), void *mctx)
{
  TaoFunctionBegin;
  if (mymonitor) {
    if (gaapp->numbermonitors >= MAX_TAO_MONITORS) {
      SETERRQ(1,"Too many monitors set.");
    }
    gaapp->monitor[gaapp->numbermonitors] = mymonitor;
    gaapp->monitorcontext[gaapp->numbermonitors++] = (void*)mctx;
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoGAMonitor"
/*@
   TaoGAAppMonitor - Apply the monitor functions for a TAO_GA_APPLICATION object

   Collective on TAO_GA_APPLICATION

   Input Parameters:
.  gaapp - the TAO_GA_APPLICATION application context

   Level: developer

.keywords: GAApplication, monitor, View

.seealso: TaoGAAppSetMonitor()
@*/
int TaoGAAppMonitor(TAO_GA_APPLICATION gaapp)
{
  int i, info;
  TaoFunctionBegin;
  for (i=0; i<gaapp->numbermonitors;i++) {
    info = (*gaapp->monitor[i])(gaapp, gaapp->monitorcontext[i]); CHKERRQ(info);

  }
  TaoFunctionReturn(0);
}



