#include "tao_general.h"
#include "taoapp_petsc.h"  /*I  "tao.h"  I*/
#include "tao.h"

static int Tao_Solve=0;
extern int TAO_APP_COOKIE;

int TaoAppGetTaoPetscApp(TAO_APPLICATION taoapp,TaoPetscApplication**tpapp);

#undef __FUNCT__  
#define __FUNCT__ "TaoSolveApplication"
/*@
  TaoSolveApplication - Find a solution to the application using the set TAO solver.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
-  tao - the TAO_SOLVER context

   Level: beginner

.keywords: solve

.seealso: TaoAppGetSolutionVec()

@*/
int TaoSolveApplication(TAO_APPLICATION taoapp, TAO_SOLVER tao){
  int info;

  PetscFunctionBegin;
  info = TaoSetupApplicationSolver(taoapp,tao);CHKERRQ(info);
  info = PetscLogEventBegin(Tao_Solve,tao,0,0,0);CHKERRQ(info);
  info = TaoSolve(tao); CHKERRQ(info);
  info = PetscLogEventEnd(Tao_Solve,tao,0,0,0);CHKERRQ(info);
  //  info = TaoGetSolutionStatus(tao,&taoapp->iter,&taoapp->fval,&taoapp->residual,0,0,&taoapp->reason); 
  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetupApplicationSolver"
/*@
   TaoSetupApplicationSolver - This routine creates the vectors,
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

   Note: 
   This routine should be called before TaoGetKSP(), but after 
   TaoAppSetInitialSolutionVec() and after TaoAppSetHessianMat() (when Newton solvers are used). 

   Note: 
   This method is called during  TaoSetOptions() and TaoSolveApplication()
   
.keywords: application, context

@*/
int TaoSetupApplicationSolver(TAO_APPLICATION myapp, TAO_SOLVER tao ){
  int info;
  TaoPetscApplication* taopetscapp;
  PetscFunctionBegin; 
  PetscValidHeaderSpecific(tao,TAO_COOKIE,2);
  if (Tao_Solve==0){
    info = PetscLogEventRegister("TaoSolve",TAO_APP_COOKIE,&Tao_Solve); CHKERRQ(info);
  }
  info = TaoAppGetTaoPetscApp(myapp,&taopetscapp);
  info = TaoSetApplication(tao,taopetscapp);CHKERRQ(info);
  taopetscapp->tao=tao;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions"
/*@
  TaoSetOptions - Sets various TAO parameters from user options

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO Application (optional)
-  tao - the TAO optimization solver (optional)
   Level: beginner

   Note: 
   This routine should be called after TaoSetupApplicationSolver()

   Note: 
   This routine must be called if there are multiple processors involved and 
   the MPI Communicator is different than MPI_COMM_WORLD.

.keywords:  options

.seealso: TaoSolveApplication()

@*/
int TaoSetOptions(TAO_APPLICATION taoapp, TAO_SOLVER tao){
  int info;
  const char *prefix=0;
  PetscTruth flg;
  MPI_Comm comm=MPI_COMM_WORLD;

  PetscFunctionBegin;

  if (tao){
    PetscValidHeaderSpecific(tao,TAO_COOKIE,2);
    info = PetscObjectGetOptionsPrefix((PetscObject)tao,&prefix); CHKERRQ(info);
    info = PetscObjectGetComm((PetscObject)tao,&comm);CHKERRQ(info);
    info = PetscOptionsBegin(comm,prefix,"TAO PETSC APPLICATIONS ","solver");CHKERRQ(info);

    info = TaoSetFromOptions(tao); CHKERRQ(info);

    flg=PETSC_FALSE;
    info = PetscOptionsName("-tao_xmonitor","Use graphics convergence","TaoPetscXMonitor",&flg);CHKERRQ(info);
    if (flg){
      info = TaoSetPetscXMonitor(tao); CHKERRQ(info);
    }

    info = PetscOptionsEnd();CHKERRQ(info);
  }

  if (taoapp){
    info = TaoAppSetFromOptions(taoapp); CHKERRQ(info);
  }

  if (tao && taoapp){
    info = TaoSetupApplicationSolver(taoapp,tao);CHKERRQ(info);
    info = PetscOptionsName("-tao_lmvmh","User supplies approximate hessian for LMVM solvers","TaoLMVMSetH0",&flg);
    if (flg){
      info=TaoBLMVMSetH0(tao,TAO_TRUE);CHKERRQ(info);
      info=TaoLMVMSetH0(tao,TAO_TRUE);CHKERRQ(info);
    }
  }
  
  PetscFunctionReturn(0);
}
 


static char TaoPetscAppXXX[] = "TaoPetscApp";

#undef __FUNCT__  
#define __FUNCT__ "TaoAppDestroyPetscAppXXX"
static int TaoAppDestroyPetscAppXXX(void*ctx){
  TaoPetscApplication *mctx=(TaoPetscApplication*)ctx;
  PetscFunctionBegin;
  if (mctx){
    delete mctx;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetTaoPetscApp"
int TaoAppGetTaoPetscApp(TAO_APPLICATION taoapp,TaoPetscApplication**tpapp){
  int info;
  PetscInt ii;
  MPI_Comm comm;
  TaoPetscApplication*ttpapp;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  info=TaoAppQueryForObject(taoapp,TaoPetscAppXXX,(void**)&ttpapp); CHKERRQ(info);
  if (ttpapp==0){
    info=PetscObjectGetComm((PetscObject)taoapp,&comm); CHKERRQ(info);
    ttpapp=new TaoPetscApplication(comm);
    info=TaoAppAddObject(taoapp,TaoPetscAppXXX,(void*)ttpapp,&ii); CHKERRQ(info);
    info=TaoAppSetDestroyRoutine(taoapp,TaoAppDestroyPetscAppXXX, (void*)ttpapp); CHKERRQ(info);
  }
  ttpapp->papp=taoapp;
  *tpapp=ttpapp;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetKSP"
/*@C
  TaoGetKSP - Gets the linear solver used by the optimization solver.
  Application writers should use TaoAppGetKSP if they need direct access
  to the PETSc KSP object.
  
   Input Parameters:
.  tao - the TAO solver

   Output Parameters:
.  ksp - the KSP linear solver used in the optimization solver

   Level: developer

.keywords: Application

.seealso: TaoAppGetKSP()

@*/
int TaoGetKSP(TAO_SOLVER tao, KSP *ksp)
{
  int info;
  TaoLinearSolver *tsolver;
  
  PetscFunctionBegin; 
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (ksp){
    *ksp=0;
    info = TaoGetLinearSolver(tao,&tsolver); CHKERRQ(info);
    if (tsolver){
      TaoLinearSolverPetsc *psolver = dynamic_cast <TaoLinearSolverPetsc *> (tsolver);
      *ksp=psolver->GetKSP();
    }
  } 
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetTaoSolver"
int TaoAppGetTaoSolver(TAO_APPLICATION taoapp, TAO_SOLVER *tao){
  int info;
  TaoPetscApplication* taopetscapp;

  PetscFunctionBegin; 
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  info = TaoAppGetTaoPetscApp(taoapp,&taopetscapp); CHKERRQ(info);
  if (tao){ *tao=taopetscapp->tao; }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetVariableBoundVecs"
/*@C
  TaoGetVariableBoundVecs - Get the vectors with the
  lower and upper bounds in current solver.

   Input Parameters:
.  tao - the TAO solver

   Output Parameter:
+  XL - the lower bounds
-  XU - the upper bounds

   Level: intermediate

   Note: 
   These vectors should not be destroyed.

.keywords: Application, bounds

.seealso: TaoAppGetGradientVec(), TaoAppGetSolutionVec(), TaoAppSetVariableBoundsRoutine()
@*/
int TaoGetVariableBoundVecs(TAO_SOLVER tao, Vec *XL, Vec *XU){
  int info;
  TaoVec *XXLL,*XXUU;
  PetscFunctionBegin; 
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = TaoGetVariableBounds(tao,&XXLL,&XXUU); CHKERRQ(info);
  if (XL){
    *XL=0;
    if (XXLL){
      TaoVecPetsc *pl = dynamic_cast <TaoVecPetsc *> (XXLL);
      *XL = pl->GetVec();
    }
  }
  if (XU){
    *XU=0;
    if (XXUU){
      TaoVecPetsc *pu = dynamic_cast <TaoVecPetsc *> (XXUU);
      *XU = pu->GetVec();
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetGradientVec"
/*@C
  TaoAppGetGradientVec - Get the vector with the
  gradient in the current application.

   Input Parameters:
.  tao - the solver

   Output Parameter:
.  G - the gradient vector

   Level: intermediate

   Note: 
   This vector should not be destroyed.

.keywords: Application, gradient

.seealso:  TaoAppGetSolutionVec()
@*/
int TaoAppGetGradientVec(TAO_SOLVER tao, Vec *G){
  int info;
  TaoVec* GG;
  PetscFunctionBegin; 
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = TaoGetGradient(tao,&GG); CHKERRQ(info);
  if (G&&GG) {
    TaoVecPetsc *pg = dynamic_cast <TaoVecPetsc *> (GG);
    *G = pg->GetVec();
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetHessianMat"
/*@C
  TaoAppGetGessianMat - Get the vector with the
  gradient in the current application.

   Input Parameters:
.  tao - the solver

   Output Parameter:
.  H - the hessian matrix
.  Hpre - the hessian preconditioner

   Level: intermediate

   Note: 
   These matrices should not be destroyed.

.keywords: Application, hessian

.seealso:  TaoAppGetSolutionVec(), TaoAppGetGradientVec()
@*/
int TaoAppGetHessianMat(TAO_SOLVER tao, Mat *H, Mat *Hpre){
  int info;
  TaoMat* TM=0;
  MatStructure flag;
  PetscFunctionBegin; 
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = TaoGetHessian(tao,&TM); CHKERRQ(info);
  if (H && TM) {
    TaoMatPetsc *tmp = dynamic_cast <TaoMatPetsc *> (TM);
    info = tmp->GetMatrix(H, Hpre, &flag); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoCopyDualsOfVariableBounds"
/*@
  TaoCopyDualsOfVariableBounds - Copy the current dual variables
  corresponding the lower and upper bounds on the variables.
  
   Input Parameters:
.  tao - the solver

   Output Parameter:
+  DXL - the lower bounds
-  DXU - the upper bounds

   Level: intermediate

   Note:  
   Existing vectors of commensurate distribution to the
   variable bounds should be passed into this routine.

   Note: 
   These variables may have to be computed.  It may not be efficient
   to call this routine in a Monitor.

   Note: 
   These variables can be interpreted as the sensitivity of
   the objective value with respect to the bounds on the variables.

.keywords: Application, bounds, duals

.seealso: TaoAppGetGradientVec(), TaoAppGetSolutionVec(), TaoAppSetVariableBoundsRoutine()
@*/
int TaoCopyDualsOfVariableBounds(TAO_SOLVER tao, Vec DXL, Vec DXU){
  int info;
  TaoVecPetsc *ddxl,*ddxu;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = TaoWrapPetscVec(DXL,&ddxl); CHKERRQ(info);
  info = TaoWrapPetscVec(DXU,&ddxu); CHKERRQ(info);
  info = TaoGetDualVariables(tao, ddxl, ddxu); CHKERRQ(info);
  info = TaoVecDestroy(ddxl); CHKERRQ(info);
  info = TaoVecDestroy(ddxu); CHKERRQ(info);
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoSetInequalityConstraints"
/*@C
   TaoSetInequality Constraints - Set inequality constraints for OOQP

   Collective on TAO_SOLVER

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  ll - vector to store lower bounds
.  uu - vector to store upper bounds
-  AA - matrix to store linear inequality constraints

   Level: intermediate

.keywords: TAO_SOLVER, inequalities

@*/
int TaoSetInequalityConstraints(TAO_APPLICATION taoapp, Vec ll, Mat A, Vec uu){

  int info;
  TaoPetscApplication* taopetscapp;
  PetscFunctionBegin;
  info = TaoAppGetTaoPetscApp(taoapp,&taopetscapp);
  info = TaoVecDestroy(taopetscapp->taofvll); CHKERRQ(info); taopetscapp->taofvll=0;
  info = TaoWrapPetscVec(ll,&taopetscapp->taofvll); CHKERRQ(info);
  info = TaoVecDestroy(taopetscapp->taofvuu); CHKERRQ(info); taopetscapp->taofvuu=0;
  info = TaoWrapPetscVec(uu,&taopetscapp->taofvuu); CHKERRQ(info);
  info = TaoMatDestroy(taopetscapp->taoAA); CHKERRQ(info); taopetscapp->taoAA=0;
  info = TaoWrapPetscMat(A,&taopetscapp->taoAA); CHKERRQ(info);
  PetscFunctionReturn(0);
}
