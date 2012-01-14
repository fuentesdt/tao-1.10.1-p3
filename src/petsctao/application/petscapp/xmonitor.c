#include "tao.h"
#include "petscdraw.h"

typedef struct{
  PetscLogDouble t0;
  PetscDraw Draw1,Draw2;
  PetscDrawLG LG1,LG2;
  TAO_SOLVER tao;
}TaoXMonitor;

int TaoPetscXMonitor(TAO_SOLVER, void *);
int DestroyTaoPetscXMonitor(void *);


int TaoSetPetscXMonitor(TAO_SOLVER tao){
  int info;
  TaoXMonitor *mctx;
  PetscFunctionBegin; 
  PetscNew(TaoXMonitor,&mctx);
  mctx->tao=tao;
  info=TaoSetMonitor(tao,TaoPetscXMonitor,(void*)mctx); CHKERRQ(info);
  info=TaoSetDestroyRoutine(tao,DestroyTaoPetscXMonitor,(void*)mctx); CHKERRQ(info);
  PetscFunctionReturn(0);
}

int TaoPetscXMonitor(TAO_SOLVER ttao, void *mctx){
  TaoXMonitor *mntr = (TaoXMonitor*) mctx;
  TAO_SOLVER  tao= mntr->tao;
  int info,rank=0;
  PetscInt iter;
  char display[51];
  PetscLogDouble t1;
  PetscReal t2,fff;
  double ff,Norm1,norm2,norm3;
  TaoTerminateReason reason;
  PetscTruth flag;

  PetscFunctionBegin; 

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  if (rank){
    PetscFunctionReturn(0);    
  }


  info = TaoGetSolutionStatus(tao,&iter,&ff,&Norm1,&norm2,&norm3,&reason);CHKERRQ(info);

  if (iter==0){
    info = PetscOptionsGetenv(PETSC_COMM_WORLD,"DISPLAY",display,50,&flag);
    CHKERRQ(info);
    info = PetscDrawCreate(PETSC_COMM_WORLD,display,
                      "Function Value",
                      0,0,500,150,&mntr->Draw1);CHKERRQ(info);
    info = PetscDrawCreate(PETSC_COMM_WORLD,display,
                      "Residual of Optimality Conditions",
                      0,175,500,150,&mntr->Draw2);CHKERRQ(info);
    info = PetscDrawSetFromOptions(mntr->Draw1); CHKERRQ(info);
    info = PetscDrawSetFromOptions(mntr->Draw2); CHKERRQ(info);
    info = PetscDrawLGCreate(mntr->Draw1,1,&mntr->LG1); CHKERRQ(info);
    info = PetscDrawLGCreate(mntr->Draw2,1,&mntr->LG2); CHKERRQ(info);
    info = PetscGetCPUTime(&mntr->t0); CHKERRQ(info);
  }

  info = PetscGetCPUTime(&t1); CHKERRQ(info);
  t1 = t1 - mntr->t0;
  t2=t1;
  fff=ff;
  info = PetscDrawLGAddPoint(mntr->LG1,&t2,&fff); CHKERRQ(info);
  info = PetscDrawLGDraw(mntr->LG1); CHKERRQ(info);
  fff=Norm1;
  info = PetscDrawLGAddPoint(mntr->LG2,&t2,&fff); CHKERRQ(info);
  info = PetscDrawLGDraw(mntr->LG2); CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DestroyTaoPetscXMonitor"
int DestroyTaoPetscXMonitor(void *ctx){

  TaoXMonitor * mntr = (TaoXMonitor *) ctx;
  int info;

  PetscFunctionBegin; 
  info = PetscSleep(-1); CHKERRQ(info);
  info =  PetscDrawLGDestroy(mntr->LG1); CHKERRQ(info);
  info =  PetscDrawLGDestroy(mntr->LG2); CHKERRQ(info);
  info =  PetscDrawDestroy(mntr->Draw1); CHKERRQ(info);
  info =  PetscDrawDestroy(mntr->Draw2); CHKERRQ(info);
  info = PetscFree(ctx); CHKERRQ(info);
  PetscFunctionReturn(0);
}
