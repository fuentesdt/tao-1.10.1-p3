#include "tao.h"
#include "petsc.h"
#include "daapp_impl.h"
#include "taodaapplication.h"

typedef struct {
  PetscLogDouble tstagebegin[20],tstageend[20];
  PetscInt nfeval[20],ngeval[20],nheval[20],nlsolves[20];
} DAAppTimeMonitor;



#undef __FUNCT__
#define __FUNCT__ "DAAppTimeMonitorBefore"
int DAAppTimeMonitorBefore(TAO_APPLICATION daapplication, DA da, PetscInt level, void *ctx){
  int info;
  DAAppTimeMonitor *dactx=(DAAppTimeMonitor*)ctx;
  PetscLogDouble t2;
  PetscFunctionBegin;
  info=PetscGetTime(&t2);CHKERRQ(info);
  dactx->tstagebegin[level]=t2;
  info=TaoAppResetCounters(daapplication);CHKERRQ(info);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "DAAppTimeMonitorAfter"
int DAAppTimeMonitorAfter(TAO_APPLICATION daapplication, DA da, PetscInt level, void *ctx){
  int info;
  PetscInt stats[4];
  DA_APPLICATION daapp;
  DAAppTimeMonitor *dactx=(DAAppTimeMonitor*)ctx;
  PetscLogDouble t0,t1,t2;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  info=PetscGetTime(&t2);CHKERRQ(info);
  dactx->tstageend[level]=t2;
  t0=t2 - dactx->tstagebegin[0];
  t1=t2 - dactx->tstagebegin[level];
  info=PetscPrintf(PETSC_COMM_WORLD," Grid %d: Time:  %4.4e \n",level,t1);
  info=PetscPrintf(PETSC_COMM_WORLD," Total Time:  %4.4e \n",t0);
  info=TaoAppCounters(daapplication,stats);CHKERRQ(info);
  dactx->nfeval[level]=stats[0];
  dactx->ngeval[level]=stats[1];
  dactx->nheval[level]=stats[2];
  dactx->nlsolves[level]=stats[3];
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppTimeMonitorAfterAll"
int DAAppTimeMonitorAfterAll(TAO_APPLICATION daapplication, DA da, PetscInt level, void *ctx){
  int info,size;
  PetscInt i,nlevels;
  PetscInt mx,my;
  DAAppTimeMonitor *dactx=(DAAppTimeMonitor*)ctx;
  PetscLogDouble t0;
  PetscFunctionBegin;
  info=DAAppGetNumberOfDAGrids(daapplication,&nlevels);CHKERRQ(info);
  if (nlevels==level+1){
    MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(info);
    PetscPrintf(PETSC_COMM_WORLD,"---------------------- \n");
    PetscPrintf(PETSC_COMM_WORLD,"Processors: %d \n",(PetscInt)size);

    PetscPrintf(PETSC_COMM_WORLD,"Mesh:                 ");
    for (i=0;i<nlevels;i++){
      info = DAAppGetDA(daapplication,i,&da);CHKERRQ(info);
      info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,
		       PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		       PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
      PetscPrintf(PETSC_COMM_WORLD," &%4d,%4d",mx,my);
    }
    PetscPrintf(PETSC_COMM_WORLD," \n");
    PetscPrintf(PETSC_COMM_WORLD,"Times:                ");
    for (i=0;i<nlevels;i++){
      t0=dactx->tstageend[i] - dactx->tstagebegin[i];
      PetscPrintf(PETSC_COMM_WORLD," & %7.2f ",t0);
    }
    PetscPrintf(PETSC_COMM_WORLD," \n");
    PetscPrintf(PETSC_COMM_WORLD,"Function Evaluations: ");
    for (i=0;i<nlevels;i++){
      PetscPrintf(PETSC_COMM_WORLD," & %7d ",dactx->nfeval[i]);
    }
    PetscPrintf(PETSC_COMM_WORLD," \n");
    PetscPrintf(PETSC_COMM_WORLD,"Gradient Evaluations: ");
    for (i=0;i<nlevels;i++){
      PetscPrintf(PETSC_COMM_WORLD," & %7d ",dactx->ngeval[i]);
    }
    PetscPrintf(PETSC_COMM_WORLD," \n");
    PetscPrintf(PETSC_COMM_WORLD,"Hessian Evaluations:  ");
    for (i=0;i<nlevels;i++){
      PetscPrintf(PETSC_COMM_WORLD," & %7d ",dactx->nheval[i]);
    }
    PetscPrintf(PETSC_COMM_WORLD," \n");
    /*
    PetscPrintf(PETSC_COMM_WORLD,"Linear Solves:        %6d ",size);
    for (i=0;i<nlevels;i++){
      PetscPrintf(PETSC_COMM_WORLD," & %4d ",dactx->nlsolves[i]);
    }
    PetscPrintf(PETSC_COMM_WORLD," \n");
    */
    PetscPrintf(PETSC_COMM_WORLD,"---------------------- \n");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppTimeMonitorDestroy"
int DAAppTimeMonitorDestroy(void *ctx){
  int info;
  PetscFunctionBegin;
  info=TaoApplicationFreeMemory(ctx);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppPrintStageTimes"
int DAAppPrintStageTimes(TAO_APPLICATION taoapp){
  int info;
  DAAppTimeMonitor *dactx;
  PetscFunctionBegin;
  PetscNew(DAAppTimeMonitor,&dactx);
  info = DAAppSetBeforeMonitor(taoapp,DAAppTimeMonitorBefore,(void*)dactx);CHKERRQ(info);
  info = DAAppSetAfterMonitor(taoapp,DAAppTimeMonitorAfter,(void*)dactx);CHKERRQ(info);
  info = DAAppSetAfterMonitor(taoapp,DAAppTimeMonitorAfterAll,(void*)dactx);CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(taoapp,DAAppTimeMonitorDestroy,(void*)dactx); CHKERRQ(info);
  PetscFunctionReturn(0);
}
  
#undef __FUNCT__
#define __FUNCT__ "DAAppInterpolationMonitor"
int DAAppInterpolationMonitor(TAO_APPLICATION taoapp, DA da, PetscInt level, void *ctx){
  
  int info;
  PetscInt mx,my,n;
  PetscScalar dd=-1.0;
  Vec X,XCoarse,XX;
  Mat Interpolate;

  PetscFunctionBegin;
  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);
  if (level>0){
    info = DAAppGetInterpolationMatrix(taoapp,level,&Interpolate,0); CHKERRQ(info);
    info = DAAppGetSolution(taoapp,level,&X); CHKERRQ(info);
    info = DAAppGetSolution(taoapp,level-1,&XCoarse); CHKERRQ(info);
    info = VecDuplicate(X,&XX); CHKERRQ(info);
    info = MatMult(Interpolate,XCoarse,XX); CHKERRQ(info);
    info = VecAXPY(XX, dd, X); CHKERRQ(info);
    info = VecNorm(XX,NORM_INFINITY,&dd); CHKERRQ(info);
    PetscPrintf(MPI_COMM_WORLD,"Maximum Interpolation Error: %4.2e\n",dd);
    info = VecNorm(XX,NORM_1,&dd); CHKERRQ(info);
    info = VecGetSize(XX,&n); CHKERRQ(info);
    PetscPrintf(MPI_COMM_WORLD,"Average Interpolation Error: %4.2e\n",dd/n);
    info = VecDestroy(XX); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "DAAppPrintInterpolationError"
int DAAppPrintInterpolationError(TAO_APPLICATION taoapp){
  int info;
  PetscFunctionBegin;
  info = DAAppSetAfterMonitor(taoapp,DAAppInterpolationMonitor,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}
