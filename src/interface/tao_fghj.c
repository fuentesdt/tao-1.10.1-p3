#include "src/tao_impl.h"      /*I "tao_solver.h"  I*/


/* --------------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoComputeFunction"
/*@C
   TaoComputeFunction - Computes the function that has been
   set with TaoSetFunction().

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  xx - input vector

   Output Parameter:
.  y - function value

   TaoComputeFunction() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, minimization, function

.seealso: TaoApplication::EvaluateObjectiveFunction(),
          TaoComputeFunctionGradient(), TaoComputeHessian()
@*/
int TaoComputeFunction(TAO_SOLVER solver,TaoVec* xx,double *y)
{
  int    info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  if ( solver->taoappl){
    info = solver->taoappl->EvaluateObjectiveFunction(xx,y); CHKERRQ(info);
  }
  solver->nfuncs++;
  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoIncrementGradientsCounter"
/*@
   TaoIncrementGradientsCounter - Increments the gradient
   counted by TAO.

   Not Collective

   Input Parameter:
+  solver -  TAO_SOLVER context
-  nevals - number of gradient evaluations to be added

   Notes:
   This counter is reset to zero for each successive call to TaoSolve().

   Level: developer

.keywords: Linear Solver, Objective Function
@*/
int TaoIncrementGradientsCounter(TAO_SOLVER solver,TaoInt nevals)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  solver->ngrads += nevals ;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoComputeGradient"
/*@C
   TaoComputeGradient - Computes the function value and its gradient

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  xx - input vector

   Output Parameter:
.  gg - gradient vector

   Notes:
   TaoComputeGradient() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, gradient

.seealso:  TaoApplication::EvaluateGradient(), TaoComputeFunction()
@*/
int TaoComputeGradient(TAO_SOLVER solver,TaoVec* xx, TaoVec* gg)
{
  int    info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  if ( solver->taoappl){
    info = solver->taoappl->EvaluateGradient(xx,gg); CHKERRQ(info);
  }
  solver->ngrads++;

  if (solver->viewgradient) {info=gg->View();CHKERRQ(info);}

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoComputeFunctionGradient"
/*@C
   TaoComputeFunctionGradient - Computes the function value and its gradient

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  xx - input vector

   Output Parameter:
+  f - function value   
-  gg - gradient vector

   Notes:
   TaoComputeFunctionGradient() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, gradient

.seealso:  TaoApplication::EvaluateObjectiveAndGradient(), TaoComputeFunction()
@*/
int TaoComputeFunctionGradient(TAO_SOLVER solver,TaoVec* xx, double *f, TaoVec* gg)
{
  int    info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  if ( solver->taoappl){
    info = solver->taoappl->EvaluateObjectiveAndGradient(xx,f,gg); CHKERRQ(info);
  }
  solver->nfgrads++;

  if (solver->viewgradient) {info=gg->View();CHKERRQ(info);}

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoComputeHessian"
/*@C
   TaoComputeHessian - Computes the Hessian matrix that has been
   set with TaoSetHessian().

   Collective on TAO_SOLVER and Mat

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  xx - input vector

   Output Parameters:
.  HH - Hessian matrix

   Notes: 
   Most users should not need to explicitly call this routine, as it
   is used internally within the minimization solvers. 

   TaoComputeHessian() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, Hessian, matrix

.seealso:  TaoApplication::EvaluateHessian(), TaoComputeFunctionGradient(),
           TaoComputeFunction()
@*/
int TaoComputeHessian(TAO_SOLVER solver,TaoVec *xx,TaoMat *HH)
{
  int    info;
  TaoTruth flg;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  if ( solver->taoappl){
    info = solver->taoappl->EvaluateHessian(xx,HH); CHKERRQ(info);
  }
  solver->nhesss++;

  info = HH->Compatible(xx,xx,&flg); CHKERRQ(info);
  if (flg==TAO_FALSE){
    SETERRQ(1,"Hessian Matrix not Compatible with solution vector");
  }

  if (solver->viewhessian){
    info=HH->View();CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoGetHessian"
/*@C
   TaoGetHessian - Sets the pointer to a TaoMat equal to the address
   a the TaoMat containing the Hessian matrix.

   Input Parameter:
+  solver - the TAO_SOLVER solver context
-  HH - address of a pointer to a TaoMat

   Output Parameters:
.  HH - address of pointer to the Hessian matrix (or TAO_NULL)

   Note:  This routine does not create a matrix.  It sets a pointer
   to the location of an existing matrix.

   Level: developer

.seealso: TaoApplication::EvaluateHessian(), TaoComputeHessian()

.keywords: TAO_SOLVER, get, Hessian
@*/
int TaoGetHessian(TAO_SOLVER solver,TaoMat **HH)
{
  int info;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  if (!solver->hessian){
    info = solver->taoappl->GetHessianMatrix(&solver->hessian); CHKERRQ(info);
  }
  if (HH)   *HH = solver->hessian;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoComputeJacobian"
/*@C
   TaoComputeJacobian - Computes the Jacobian matrix that has been
   set.

   Collective on TAO_SOLVER and Mat

   Input Parameters:
+  solver - the TAO_SOLVER solver context
.  xx - input vector

   Output Parameters:
+  JJ - Jacobian matrix

   Notes: 
   TaoComputeJacobian() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, Jacobian, matrix

.seealso:  TaoApplication::EvaluateJacobian()
@*/
int TaoComputeJacobian(TAO_SOLVER solver,TaoVec* xx,TaoMat *JJ)
{
  int    info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  if ( solver->taoappl){
    info = solver->taoappl->EvaluateJacobian(xx,JJ); CHKERRQ(info);
  }
  solver->njac++;

  if (solver->viewjacobian){ info=JJ->View();CHKERRQ(info);}

  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoEvaluateVariableBounds"
/*@
   TaoEvaluateVariableBounds - Evaluate the  lower and upper bounds on the variables.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
.  xxll - vector of lower bounds upon the solution vector
-  xxuu - vector of upper bounds upon the solution vector

   Level: developer

.keywords: bounds

.seealso: TaoGetVariableBounds()
@*/
int TaoEvaluateVariableBounds(TAO_SOLVER tao,TaoVec *xxll,TaoVec *xxuu)
{
  int info;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = tao->taoappl->EvaluateVariableBounds(xxll,xxuu); CHKERRQ(info);
  info = TaoCheckBounds(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}


/*
#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetConstraintsBounds"
int TaoAppSetConstraintsBounds(TAO_SOLVER solver,TaoVec *RXL,TaoMat *AA, TaoVec *RXU)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  solver->RXL=RXL;
  solver->CA=AA;
  solver->RXU=RXU;
  TaoFunctionReturn(0);
}
*/

#undef __FUNCT__  
#define __FUNCT__ "TaoComputeConstraints"
/*@C
   TaoComputeConstraints - Computes the constraint vector that has been
   set in the application.

   Collective on TAO_SOLVER and Mat

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  xx - input vector

   Output Parameters:
.  rr - Constraint values

   Notes: 
   TaoComputeConstraints() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, constraint

.seealso:  TaoAppSetConstraintRoutine()
@*/
int TaoComputeConstraints(TAO_SOLVER solver,TaoVec* xx,TaoVec* rr)
{
  int    info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  if ( solver->taoappl){
    info = solver->taoappl->EvaluateConstraints(xx,rr); CHKERRQ(info);
  }
  solver->nvfunc++;
  if (solver->viewvfunc){ info=rr->View(); CHKERRQ(info); }

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetConstraints"
/*@C
   TaoGetConstraints - Get the constraint vector that has been
   set.

   Collective on TAO_SOLVER and Mat

   Input Parameters:
.  solver - the TAO_SOLVER solver context

   Output Parameters:
.  rr - Constraint values

   Level: developer

.keywords: TAO_SOLVER, compute, constraints

.seealso:  TaoApplication::EvaluateConstraints()
@*/
int TaoGetConstraints(TAO_SOLVER solver, TaoVec** rr)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  *rr = solver->vfunc;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetJacobian"
/*@C
   TaoGetJacobian - Get the Jacobian matrix

   Collective on TAO_SOLVER and Mat

   Input Parameters:
.  solver - the TAO_SOLVER solver context

   Output Parameters:
.  JJ - Jacobian Matrix

   Level: advanced

.keywords: TAO_SOLVER, compute, constraint

.seealso:  TaoApplication::EvaluateJacobian(), TaoAppSetConstraintRoutine()
@*/
int TaoGetJacobian(TAO_SOLVER solver, TaoMat **JJ)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  *JJ = solver->jacobian;
  TaoFunctionReturn(0);
}




#undef __FUNCT__  
#define __FUNCT__ "TaoGetSolution"
/*@C
   TaoGetSolution - Sets a pointer to a TaoVec to the address of the 
   vector containing the current solution.

   Input Parameter:
+  solver - the TAO_SOLVER solver context
-  xx - the address of a pointer to a TaoVec

   Output Parameter:
.  xx - the solution

   Level: advanced
   
   Note:  This routine does not create a vector.  It sets a pointer
   to the location of an existing vector.

   Note:
   This vector is a reference to the vector set in the application and passed
   to TAO in TaoSetApplication().

.keywords: solve, solution

.seealso: TaoCreate(), TaoGetGradient(), TaoSetApplication()
@*/
int TaoGetSolution(TAO_SOLVER solver,TaoVec** xx)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  *xx=solver->vec_sol;
  TaoFunctionReturn(0);
}  


#undef __FUNCT__  
#define __FUNCT__ "TaoSetLagrangianGradientVector"
/*@
   TaoSetLagrangianGradientVector - Sets a pointer to the address of a TaoVector that
   contains the gradient of the Lagrangian function.

   Input Parameter:
+  solver - the TAO_SOLVER solver context
-  gg - the gradient of the Lagrangian function

   Level: developer

   Note:  This routine does not create a vector.  The vector specified
   here will be returned whenever TaoGetGradient() is called.

.keywords: Gradient

.seealso: TaoGetGradient()

@*/
int TaoSetLagrangianGradientVector(TAO_SOLVER solver,TaoVec* gg)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  solver->vec_grad =gg;
  TaoFunctionReturn(0);
} 

#undef __FUNCT__  
#define __FUNCT__ "TaoGetGradient"
/*@C
   TaoGetGradient - Sets a pointer to the address of a TaoVector that
   contains the gradient of the Lagrangian function.

   Input Parameter:
.  solver - the TAO_SOLVER solver context

   Output Parameter:
.  gg - the gradient of the Lagrangian function

   Level: advanced

   Note:  This routine does not create a vector.  It sets a pointer
   to the location of an existing vector.

.keywords: Gradient

.seealso: TaoGetSolution(), TaoGetSolutionStatus(), TaoGetHessian(), TaoSetApplication()

@*/
int TaoGetGradient(TAO_SOLVER solver,TaoVec** gg)
{
  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  *gg = solver->vec_grad;
  TaoFunctionReturn(0);
} 


/* --------------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoComputeMeritFunction"
/*@C
   TaoComputeMeritFunction - Computes the function that has been
   set with TaoSetMeritFunction().

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  xx - input vector

   Output Parameter:
.  y - merit function value

   TaoComputeMeritFunction() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, minimization, function

.seealso: TaoSetMeritFunction()
@*/
int TaoComputeMeritFunction(TAO_SOLVER solver,TaoVec* xx,double *y)
{
  int    info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  if ( solver->MeritFunctionApply){
    info = solver->MeritFunctionApply(solver,xx,y,solver->meritctx); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetMeritFunction"
/*@C
   TaoSetMeritFunction - Sets the routine that evaluates the merit
   function.

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
.  func - merit function evaluation routine
.  funcgrad - merit function and gradient evalution routine
.  grad - merit gradient evaluation routine
.  destroy - context destructor routine
-  ctx - [optional] user-defined context for private data for the 
         function and gradient evaluation routine (may be TAO_NULL)

   Calling sequence of func:
$     func (TAO_SOLVER solver,TaoVec *xx,double *f,void *ctx);

+  solver - the TAO_SOLVER solver context
.  xx - variable vector
.  f - objective function value
-  ctx - [optional] user-defined function context 

   Calling sequence of funcgrad:
$     funcgrad (TAO_SOLVER solver,TaoVec *xx,double *f,TaoVec *gg, void *ctx);

+  solver - the TAO_SOLVER solver context
.  xx - variable vector
.  f - objective function value
.  gg - gradient vector
-  ctx - [optional] user-defined function context 

   Calling sequence of grad:
$     grad (TAO_SOLVER solver,TaoVec *xx,TaoVec *gg,void *ctx);

+  solver - the TAO_SOLVER solver context
.  xx - variable vector
.  gg - gradient vector
-  ctx - [optional] user-defined function context 

   Calling sequence of destroy:
$     destroy (TAO_SOLVER solver, void *ctx);

+  solver - the TAO_SOLVER solver context
-  ctx - [optional] user-defined function context 

   Level: developer

.keywords: TAO_SOLVER, merit function

.seealso:  TaoComputeMeritFunction()
@*/
int TaoSetMeritFunction(TAO_SOLVER solver,int (*func)(TAO_SOLVER,TaoVec*,double*,void*),
			int (*funcgrad)(TAO_SOLVER,TaoVec*,double*,TaoVec*,void*),
			int (*grad)(TAO_SOLVER,TaoVec*,TaoVec*,void*),
			int (*funcgts)(TAO_SOLVER,TaoVec*,TaoVec*,double*,double*,void*),
			int (*destroy)(TAO_SOLVER,void*),
			void *ctx)
{
  int info;
  TaoFunctionBegin;

  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);
  info = TaoMeritFunctionDestroy(solver); CHKERRQ(info);
  solver->MeritFunctionApply =func;
  solver->MeritFunctionGradientApply =funcgrad;
  solver->MeritGradientApply =grad;
  solver->MeritFunctionGTSApply = funcgts;
  solver->MeritFunctionDestroy =destroy;
  solver->meritctx           = ctx;

  TaoFunctionReturn(0);
}

/* --------------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoComputeMeritFunctionGradient"
/*@C
   TaoComputeMeritFunctionGradient - Computes the function that has been
   set with TaoSetMeritFunction().

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  xx - input vector

   Output Parameter:
.  y - merit function value
.  gg - merit gradient vector

   TaoComputeMeritFunctionGradient() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, minimization, function

.seealso: TaoSetMeritFunctionGradient()
@*/
int TaoComputeMeritFunctionGradient(TAO_SOLVER solver,TaoVec* xx,double *y,TaoVec*gg)
{
  int    info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);

  if ( solver->MeritFunctionGradientApply){
    info = solver->MeritFunctionGradientApply(solver,xx,y,gg,solver->meritctx); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}


/* --------------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoComputeMeritGradient"
/*@C
   TaoComputeMeritGradient - Computes the function that has been
   set with TaoSetMeritFunction().

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
-  xx - input vector

   Output Parameter:
.  gg - merit gradient vector

   TaoComputeMeritFunctionGradient() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, minimization, function

.seealso: TaoSetMeritFunction()
@*/
int TaoComputeMeritGradient(TAO_SOLVER solver,TaoVec* xx,TaoVec*gg)
{
  int    info;
  double ignore;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);


  if ( solver->MeritGradientApply){
    info = solver->MeritGradientApply(solver,xx,gg,solver->meritctx); CHKERRQ(info);
  } else if (solver->MeritFunctionGradientApply) {
    info = solver->MeritFunctionGradientApply(solver,xx,&ignore,gg,solver->meritctx); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

/* --------------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoComputeMeritFunctionGTS"
/*@C
   TaoComputeMeritFunctionGTS - Computes the function that has been
   set with TaoSetMeritFunction() and the inner product of the gradient
   and the step direction.

   Collective on TAO_SOLVER

   Input Parameters:
+  solver - the TAO_SOLVER solver context
.  xx - input vector
-  ss - step direction

   Output Parameter:
+  f - the function value at xx
-  gts - the inner product of the gradient and the step diretion

   TaoComputeMeritFunctionGTS() is typically used within minimization
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.keywords: TAO_SOLVER, compute, minimization, function

.seealso: TaoSetMeritFunction()
@*/
int TaoComputeMeritFunctionGTS(TAO_SOLVER solver,TaoVec* xx,TaoVec *ss, double *f,double *gts)
{
  int    info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(solver,TAO_COOKIE,1);


  if ( solver->MeritFunctionGTSApply){
    info = solver->MeritFunctionGTSApply(solver,xx,ss,f,gts,solver->meritctx); CHKERRQ(info);
  } else {
    SETERRQ(1,"MeritFunctionGTS not set");
  }

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoMeritFunctionDestroy"
/*@C
   TaoMeritFunctionDestroy - Destroy the data structures associated with
   the merit function and set associated function pointers to NULL.

   Collective on TAO_SOLVER

   Input Parameters:
.  solver - the TAO_SOLVER solver context


   Level: developer

.keywords: TAO_SOLVER, compute, minimization, function

.seealso: TaoSetMeritFunction()
@*/
int TaoMeritFunctionDestroy(TAO_SOLVER solver){
  int info;
  TaoFunctionBegin;
  if (solver->MeritFunctionDestroy){
    info=solver->MeritFunctionDestroy(solver,solver->meritctx); CHKERRQ(info);
  }
  solver->MeritFunctionApply=0;
  solver->MeritFunctionGradientApply=0;
  solver->MeritGradientApply=0;
  solver->MeritFunctionDestroy=0;
  solver->meritctx=0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetApplication"
/*@C
   TaoSetApplication - Sets the user defined context for 
   use by the optimization solver.  The application provides
   the solver with function and derivative information as
   well as data structures uses to store this information.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER solver context
-  taoapp - user application context

   Note:
   For PETSc users the TaoApplication* object is actually
   TAO_APPLICATION structure.

   
   Level: advanced

.keywords: application, context

@*/
int TaoSetApplication(TAO_SOLVER tao, TaoApplication *myapp){
  int info;
  TaoVec *xx,*RR;
  TaoMat *HH, *JJ;
  TaoVec *RXL, *RXU;
  TaoMat *CA;

  TaoFunctionBegin;
  tao->taoappl=myapp;
  info = myapp->GetVariableVector(&xx);CHKERRQ(info);
  info = myapp->GetHessianMatrix(&HH);CHKERRQ(info);
  info = myapp->GetJacobianMatrix(&JJ); CHKERRQ(info);
  info = myapp->GetConstraintVector(&RR);CHKERRQ(info);
  info = myapp->GetInequalityConstraints(&RXL,&CA,&RXU);

  tao->hessian=HH;
  tao->vec_sol=xx;
  tao->jacobian=JJ;
  tao->vfunc=RR;
  tao->RXL=RXL;
  tao->RXU=RXU;
  tao->CA=CA;

  info = TaoSetUp(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

class TaoH0Mat: public TaoMat{
 protected:
  TaoApplication *H0;
 public:  
  TaoH0Mat(TaoApplication*);
  ~TaoH0Mat();
  int Solve(TaoVec*, TaoVec*, TaoTruth*);
};

#undef __FUNCT__
#define __FUNCT__ "TaoH0Mat::TaoH0Mat"
TaoH0Mat::TaoH0Mat(TaoApplication* theappobject){
  this->H0=theappobject;
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoH0Mat::~TaoH0Mat"
TaoH0Mat::~TaoH0Mat(){
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoH0Mat::Solve"
int TaoH0Mat::Solve(TaoVec* tb, TaoVec* dx, TaoTruth *tt){
  int info;
  info=this->H0->HessianSolve(tb,dx,tt); CHKERRQ(info);
  return 0;
}

#include "src/unconstrained/impls/lmvm/lmvm.h"
#undef __FUNCT__
#define __FUNCT__ "TaoLMVMSetH0"
int TaoLMVMSetH0(TAO_SOLVER tao, TaoTruth flag)
{
  TAO_LMVM *lmvm;
  TaoMat *H0;
  int info;

  TaoFunctionBegin;
  info = TaoGetSolverContext(tao, "tao_lmvm", (void **)&lmvm); CHKERRQ(info);
  if (lmvm && lmvm->M) {
    if (TAO_TRUE == flag) {
      H0 = new TaoH0Mat(tao->taoappl);
      info = lmvm->M->SetH0(H0); CHKERRQ(info);
    }
    else {
      info = lmvm->M->SetH0(0); CHKERRQ(info);
    }
  }
  TaoFunctionReturn(0);
}

#include "src/bound/impls/blmvm/blmvm.h"
#undef __FUNCT__
#define __FUNCT__ "TaoBLMVMSetH0"
int TaoBLMVMSetH0(TAO_SOLVER tao, TaoTruth flag)
{
  TAO_BLMVM *blmvm;
  TaoMat *H0;
  int info;

  TaoFunctionBegin;
  info = TaoGetSolverContext(tao, "tao_blmvm", (void **)&blmvm); CHKERRQ(info);
  if (blmvm && blmvm->M) {
    if (TAO_TRUE == flag) {
      H0 = new TaoH0Mat(tao->taoappl);
      info = blmvm->M->SetH0(H0); CHKERRQ(info);
    }
    else {
      info = blmvm->M->SetH0(0); CHKERRQ(info);
    }
  }
  TaoFunctionReturn(0);
}

