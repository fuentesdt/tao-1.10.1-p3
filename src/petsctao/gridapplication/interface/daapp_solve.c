#include "taodaapplication.h"   /*I  "taodaapplication.h"  I*/
//#include "taoapp.h"
#include "tao.h"
#include "petsc.h"
#include "src/petsctao/application/petscapp/tao_app_impl.h"
#include "daapp_impl.h"

static int Tao_DA_Solve=0;
extern int DAAPP_COOKIE;
int TaoAppDAApp(TAO_APPLICATION, DA_APPLICATION *);

#undef __FUNCT__
#define __FUNCT__ "TaoDAAppSolve"
/*@
  TaoDAAppSolve - Solve the PETSC DA application.

   Input Parameters:
.  daapplication - the TAO DAApplication structure

   Level: beginner

.keywords: Application, DA, Solve
@*/
int TaoDAAppSolve(TAO_APPLICATION daapplication, TAO_SOLVER tao){
  int info;
  PetscInt i,j,iter;
  PetscInt mx,my,mz;
  double fval,residual;
  TaoTerminateReason reason;
  DA_APPLICATION daapp;

  PetscFunctionBegin;

  PetscValidHeaderSpecific(tao,TAO_COOKIE,2);
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);

  if (Tao_DA_Solve==0){
    info=PetscLogEventRegister("TaoSolveDAApp",DAAPP_COOKIE,&Tao_DA_Solve); CHKERRQ(info);
  }

  info = PetscLogEventBegin(Tao_DA_Solve,tao,daapp,0,0);CHKERRQ(info);
  for (i=0;i<daapp->nda; i++){
    
    
    info=TaoAppResetCounters(daapplication);CHKERRQ(info);
    info = DAGetInfo(daapp->grid[i].da,PETSC_NULL,&mx,&my,&mz,PETSC_NULL,
		     PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		     PETSC_NULL,PETSC_NULL);CHKERRQ(info);

    info = PetscInfo2(daapp->grid[i].da,"Level %d of %d in TAO_DA_APPLICATION object.\n",i,daapp->nda); CHKERRQ(info);
    info = PetscInfo3(daapp->grid[i].da,"Global dimensions of DA:  mx=%d, my=%d, mz=%d .\n",mx,my,mz); CHKERRQ(info);

    daapp->currentlevel=i;
    PetscValidHeaderSpecific(daapp->grid[i].da,DM_COOKIE,1);
    info=TaoAppSetColoring(daapplication, daapp->grid[i].coloring);CHKERRQ(info);

    if (i>0){
      info = TaoSetDown(tao); CHKERRQ(info);
      info = MatMult(daapp->grid[i].Interpolate,daapp->grid[i-1].X,daapp->grid[i].X); CHKERRQ(info);
    } 

    info = TaoAppSetInitialSolutionVec(daapplication,daapp->grid[i].X);CHKERRQ(info);

    if (daapp->grid[i].H){
      if (!daapplication->ksp) {
	MPI_Comm comm;
	info = PetscObjectGetComm((PetscObject)daapp->grid[i].H,&comm); CHKERRQ(info);
	info = KSPCreate(comm,&daapp->grid[i].ksp); CHKERRQ(info);
	info = KSPSetFromOptions(daapp->grid[i].ksp); CHKERRQ(info);
	daapplication->ksp = daapp->grid[i].ksp;
	/*
	if (tao->ksp) {
	  tao->ksp->SetKSP(daapp->grid[i].ksp);
	}
	info = TaoWrapKSP(daapp->grid[i].ksp,&tao->ksp); CHKERRQ(info);
	*/
      }
      if (daapp->IsComplementarity==PETSC_FALSE){
	info = TaoAppSetHessianMat(daapplication,daapp->grid[i].H,daapp->grid[i].H);CHKERRQ(info);
      } else {
	info = TaoAppSetJacobianMat(daapplication,daapp->grid[i].H,daapp->grid[i].H);CHKERRQ(info);
      }
    }

    info = TaoSetupApplicationSolver(daapplication, tao); CHKERRQ(info);
    for (j=0;j<daapp->nbeforemonitors;j++){
      info = PetscInfo(daapp->grid[i].da,"TAO: Call before user monitor for DA_APPLICATION object.\n"); CHKERRQ(info);
      info = (*daapp->beforemonitor[j])(daapplication,daapp->grid[i].da,i,
					daapp->beforemonitorctx[j]); CHKERRQ(info);
    }

    //        info = TaoSetUp(tao); CHKERRQ(info);

    info = PetscInfo1(daapp->grid[i].da,"TAO: Begin solving level %d of DA_APPLICATION object.\n",i); CHKERRQ(info);
    info = TaoSolveApplication(daapplication,tao);CHKERRQ(info);

    info = TaoGetSolutionStatus(tao,&iter,&fval,&residual,0,0,&reason); CHKERRQ(info);
    if (reason <= 0 ){
      info = PetscInfo1(daapp->grid[i].da,"FAILURE TO FIND SOLUTION at level %d of DA_APPLICATION.\n",i+1); CHKERRQ(info);    
      info = PetscInfo1(daapp->grid[i].da,"  TAO Reason for termination: %d.\n",reason); CHKERRQ(info);
    } else {
      info = PetscInfo1(daapp->grid[i].da,"Found solution to DA_APPLICATION at level: %d.\n",i+1);CHKERRQ(info);
      info = PetscInfo3(daapp->grid[i].da,"Iterations: %d, Objective Value: %10.8e, Residual: %4.2e.\n",
			   iter,fval,residual);CHKERRQ(info);
    }

    for (j=0;j<daapp->naftermonitors;j++){
      info = PetscInfo(daapp->grid[i].da,"TAO: Call after user monitor for DA_APPLICATION object.\n"); CHKERRQ(info);
      info = (*daapp->aftermonitor[j])(daapplication,daapp->grid[i].da,i,
				       daapp->aftermonitorctx[j]); CHKERRQ(info);
    }

    if (daapplication->ksp) {
      info = KSPDestroy(daapplication->ksp); CHKERRQ(info);
      daapplication->ksp=0;
    }
    if (i < daapp->nda-1){
      info = TaoSetDown(tao); CHKERRQ(info);
    } 
    
  }

  info = PetscLogEventEnd(Tao_DA_Solve,tao,daapp,0,0);CHKERRQ(info);
  PetscFunctionReturn(0);
}

 
