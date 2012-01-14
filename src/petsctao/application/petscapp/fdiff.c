#include "taoapp.h"         /*I  "taoapp.h"  I*/
#include "src/tao_impl.h"      /*I "src/tao_impl.h"  I*/
#include "src/petsctao/include/taopetsc.h"          /*I "src/petsctao/include/taopetsc.h" I*/
#include "petscsnes.h"

extern int TaoAppSetColoring(TAO_APPLICATION, ISColoring);
extern int TaoAppGetColoring(TAO_APPLICATION, ISColoring*);


#undef __FUNCT__  
#define __FUNCT__ "Ftemp"
static int Ftemp(SNES snes ,Vec X,Vec G,void*ctx){
  int info;
  TAO_APPLICATION taoapp = (TAO_APPLICATION)ctx;
  TAO_SOLVER tao;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ctx,TAO_APP_COOKIE,4);
  info = TaoAppGetTaoSolver(taoapp,&tao);
  info=TaoAppComputeGradient(taoapp,X,G); CHKERRQ(info);
  tao->ngrads++;
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoAppDefaultComputeGradient"
/*@C
  TaoAppDefaultComputeHessian - computes the gradient using finite differences.
 
  Collective on TAO_APPLICATION

  Input Parameters:
+ taoapp - the TAO_APPLICATION context
. X - compute gradient at this point
- ctx - the TAO_APPLICATION structure, cast to (void*)

  Output Parameters:
. G - Gradient Vector

   Options Database Key:
+  -tao_fd_gradient - Activates TaoAppDefaultComputeGradient()
-  -tao_fd_delta <delta> - change in x used to calculate finite differences




   Level: intermediate

   Note:
   This routine is slow and expensive, and is not currently optimized
   to take advantage of sparsity in the problem.  Although
   TaoAppDefaultComputeGradient is not recommended for general use
   in large-scale applications, It can be useful in checking the
   correctness of a user-provided gradient.  Use the tao method "tao_fd_test"
   to get an indication of whether your gradient is correct.


   Note:
   The gradient evaluation must be set using the routine TaoAppSetGradientRoutine().

.keywords: TAO_APPLICATION, finite differences, Hessian

.seealso: TaoAppDefaultComputeGradient(),  TaoAppSetGradientRoutine()

@*/
int TaoAppDefaultComputeGradient(TAO_APPLICATION taoapp,Vec X,Vec G,void*) 
{
  Vec TempX;
  double *g;
  double f, f2;
  int info;
  PetscInt low,high,N,i;
  PetscTruth flg;
  double h=1.0e-6;
  TAO_SOLVER tao;
  PetscFunctionBegin;
  info = TaoAppGetTaoSolver(taoapp, &tao);
  //PetscStackPush("TAO Finite Difference gradient");
  //info = PetscLogEventBegin(Tao_FiniteDifferenceGradient,taoapp,X,0,0);
  info = TaoAppComputeObjective(taoapp, X,&f); CHKERRQ(info);
  tao->nfuncs++;
  info = PetscOptionsGetReal(PETSC_NULL,"-tao_fd_delta",&h,&flg);
  info = VecDuplicate(X,&TempX); CHKERRQ(info);
  info = VecCopy(X,TempX); CHKERRQ(info);
  info = VecGetSize(X,&N); CHKERRQ(info);
  info = VecGetOwnershipRange(TempX,&low,&high); CHKERRQ(info);
  info = VecGetArray(G,&g); CHKERRQ(info);
  for (i=0;i<N;i++) {
    info = VecSetValue(TempX,i,h,ADD_VALUES);
    info = TaoAppComputeObjective(taoapp,TempX,&f2); CHKERRQ(info);
    tao->nfuncs++;
    info = VecSetValue(TempX,i,-h,ADD_VALUES);
    if (i>=low && i<high) {
      g[i-low]=(f2-f)/h;
    }

  }
  info = VecRestoreArray(G,&g); CHKERRQ(info);
  
  //info = PetscLogEventEnd(Tao_FiniteDifferenceGradient, taoapp, X, 0,0); CHKERRQ(info);
  //PetscStackPop;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppDefaultComputeHessian"
/*@C
   TaoAppDefaultComputeHessian - Computes the Hessian using finite differences. 

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context 
.  V - compute Hessian at this point
-  ctx - the TAO_APPLICATION structure, cast to (void*)

   Output Parameters:
+  H - Hessian matrix (not altered in this routine)
.  B - newly computed Hessian matrix to use with preconditioner (generally the same as H)
-  flag - flag indicating whether the matrix sparsity structure has changed

   Options Database Key:
+  -tao_fd - Activates TaoAppDefaultComputeHessian()
-  -tao_view_hessian - view the hessian after each evaluation using PETSC_VIEWER_STDOUT_WORLD

   Level: intermediate

   Notes:
   This routine is slow and expensive, and is not currently optimized
   to take advantage of sparsity in the problem.  Although
   TaoAppDefaultComputeHessian() is not recommended for general use
   in large-scale applications, It can be useful in checking the
   correctness of a user-provided Hessian.

   Note:
   The gradient evaluation must be set using the routine TaoAppSetGradientRoutine().

.keywords: TAO_APPLICATION, finite differences, Hessian

.seealso: TaoAppSetHessianRoutine(), TaoAppDefaultComputeHessianColor(), SNESDefaultComputeJacobian(),
          TaoAppSetGradientRoutine(), TaoAppDefaultComputeGradient()

@*/
int TaoAppDefaultComputeHessian(TAO_APPLICATION taoapp,Vec V,Mat *H,Mat *B,
			     MatStructure *flag,void *ctx){
  int                  info;
  MPI_Comm             comm;
  Vec                  G;
  SNES                 snes;
  TAO_SOLVER tao;


  PetscFunctionBegin;
  info = TaoAppGetTaoSolver(taoapp, &tao);

  PetscValidHeaderSpecific(V,VEC_COOKIE,2);
  info = VecDuplicate(V,&G);CHKERRQ(info);

  info = PetscInfo(G,"TAO Using finite differences w/o coloring to compute matrix\n"); CHKERRQ(info);

  info = TaoAppComputeGradient(taoapp,V,G); CHKERRQ(info);
  tao->ngrads++;

  info = PetscObjectGetComm((PetscObject)(*H),&comm);CHKERRQ(info);
  info = SNESCreate(comm,&snes);CHKERRQ(info);

  info = SNESSetFunction(snes,G,Ftemp,taoapp);CHKERRQ(info);
  info = SNESDefaultComputeJacobian(snes,V,H,B,flag,taoapp);CHKERRQ(info);

  info = SNESDestroy(snes);CHKERRQ(info);
  
  info = VecDestroy(G);CHKERRQ(info);
  
  PetscFunctionReturn(0);
}

EXTERN_C_END



EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoAppDefaultComputeHessianColor"
/*@C
   TaoAppDefaultComputeHessianColor - Computes the Hessian using colored finite differences. 

   Collective on TAO_APPLICATION

   Input Parameters:
+  tao - the TAO_APPLICATION context
.  V - compute Hessian at this point
-  ctx - the TAO_APPLICATION structure, cast to (void*)

   Output Parameters:
+  H - Hessian matrix (not altered in this routine)
.  B - newly computed Hessian matrix to use with preconditioner (generally the same as H)
-  flag - flag indicating whether the matrix sparsity structure has changed

   Options Database Keys:
+  -mat_fd_coloring_freq <freq>
-  -tao_view_hessian - view the hessian after each evaluation using PETSC_VIEWER_STDOUT_WORLD

   Level: intermediate

   Note:
   The gradient evaluation must be set using the routine TaoSetPetscGradient().

 .keywords: TAO_APPLICATION, finite differences, Hessian, coloring, sparse

.seealso: TaoAppSetHessianRoutine(), TaoAppDefaultComputeHessian(),SNESDefaultComputeJacobianColor(), 
          TaoAppSetGradientRoutine(), TaoAppSetColoring()

@*/
int TaoAppDefaultComputeHessianColor(TAO_APPLICATION taoapp, Vec V, Mat *HH,Mat *BB,
				  MatStructure *flag,void *ctx){
  int                 info;
  MPI_Comm            comm;
  Vec                 G=0;
  Mat                 H=*HH,B=*BB;
  SNES                snes;
  ISColoring          iscoloring;
  MatFDColoring       matfdcoloring;
  TAO_SOLVER          tao;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(H,MAT_COOKIE,3);  
  PetscValidHeaderSpecific(B,MAT_COOKIE,4);  
  PetscCheckSameComm(V,2,H,3);
  PetscCheckSameComm(H,3,B,4);

  info = TaoAppGetTaoSolver(taoapp,&tao); CHKERRQ(info);
  info = TaoAppGetColoring(taoapp,&iscoloring); CHKERRQ(info);
  if (!iscoloring){
    SETERRQ(1,"Must set coloring before using this routine.  Try Finite Differences without coloring\n");
  }
  info = VecDuplicate(V,&G);CHKERRQ(info);

  info=PetscInfo(G,"TAO computing matrix using finite differences and coloring\n"); CHKERRQ(info);

  info=TaoAppComputeGradient(taoapp,V,G); CHKERRQ(info);
  tao->ngrads++;

  info = PetscObjectGetComm((PetscObject)(H),&comm);CHKERRQ(info);
  info = SNESCreate(comm,&snes);CHKERRQ(info);

  info = MatFDColoringCreate(H,iscoloring,&matfdcoloring);CHKERRQ(info);
  info = MatFDColoringSetFunction(matfdcoloring,(int (*)(void)) Ftemp,taoapp);CHKERRQ(info);
  info = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(info);

  info = SNESSetFunction(snes,G,Ftemp,taoapp);CHKERRQ(info);
  info = SNESSetJacobian(snes,H,B,SNESDefaultComputeJacobianColor,(void*)matfdcoloring);CHKERRQ(info);
  info = SNESDefaultComputeJacobianColor(snes,V,HH,BB,flag,matfdcoloring);CHKERRQ(info);

  info = MatFDColoringDestroy(matfdcoloring);CHKERRQ(info);
  info = SNESDestroy(snes);CHKERRQ(info);
  
  info = VecDestroy(G);CHKERRQ(info);
  PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetFiniteDifferencesOptions"
/* @
  TaoAppSetFiniteDifferencesOptions - Sets various TAO parameters from user options

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO Application (optional)

   Level: beginner

.keywords:  options, finite differences

.seealso: TaoSolveApplication();

@ */
int TaoAppSetFiniteDifferencesOptions(TAO_APPLICATION taoapp){
  int info;
  PetscTruth flg;

  PetscFunctionBegin;

  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);

  flg=PETSC_FALSE;
  info = PetscOptionsName("-tao_fd","use finite differences for Hessian","TaoAppDefaultComputeHessian",&flg);CHKERRQ(info);
  if (flg) {
    info = TaoAppSetHessianRoutine(taoapp,TaoAppDefaultComputeHessian,(void*)taoapp);CHKERRQ(info);
    info = PetscInfo(taoapp,"Setting default finite difference Hessian matrix\n"); CHKERRQ(info);
  }

  flg=PETSC_FALSE;
  info = PetscOptionsName("-tao_fdgrad","use finite differences for gradient","TaoAppDefaultComputeGradient",&flg);CHKERRQ(info);
  if (flg) {
    info = TaoAppSetGradientRoutine(taoapp,TaoAppDefaultComputeGradient,(void*)taoapp);CHKERRQ(info);
    info = PetscInfo(taoapp,"Setting default finite difference gradient routine\n"); CHKERRQ(info);
  }


  flg=PETSC_FALSE;
  info = PetscOptionsName("-tao_fd_coloring","use finite differences with coloring to compute Hessian","TaoAppDefaultComputeHessianColor",&flg);CHKERRQ(info);
  if (flg) {
    info = TaoAppSetHessianRoutine(taoapp,TaoAppDefaultComputeHessianColor,(void*)taoapp);CHKERRQ(info);
    info = PetscInfo(taoapp,"Use finite differencing with coloring to compute Hessian \n"); CHKERRQ(info);
  }
    
  PetscFunctionReturn(0);
}


static char TaoAppColoringXXX[] = "TaoColoring";

typedef struct {
  ISColoring coloring;
} TaoAppColorit;

#undef __FUNCT__  
#define __FUNCT__ "TaoAppDestroyColoringXXX"
static int TaoAppDestroyColoringXXX(void*ctx){
  int info;
  TaoAppColorit *mctx=(TaoAppColorit*)ctx;
  PetscFunctionBegin;
  if (mctx){
    /*
    if (mctx->coloring){  
      info = ISColoringDestroy(mctx->coloring);CHKERRQ(info);
    }
    */
    info = PetscFree(mctx); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetColoring"
/*@
   TaoAppSetColoring - Set the matrix coloring to be used when computing the
   Hessian by finite differences.

   Collective on TAO_APPLICATION

   Input Parameters:
+  app - the TAO_APPLICATION context
-  coloring - the coloring to be used

   Level: intermediate

   Note:
   The gradient evaluation must be set using the routine TaoSetPetscGradient().

 .keywords: TAO_APPLICATION, finite differences, Hessian, coloring, sparse

.seealso: TaoAppSetHessianRoutine(), TaoAppDefaultComputeHessianColor(),
          TaoAppSetGradientRoutine()

@*/
int TaoAppSetColoring(TAO_APPLICATION taoapp, ISColoring coloring){
  int info;
  PetscInt ii;
  TaoAppColorit *mctx;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  
  info = TaoAppQueryForObject(taoapp,TaoAppColoringXXX,(void**)&mctx); CHKERRQ(info);
  if (mctx==0){
    info=PetscNew(TaoAppColorit,(void**)&mctx); CHKERRQ(info);
    info=TaoAppAddObject(taoapp,TaoAppColoringXXX,(void*)mctx,&ii); CHKERRQ(info);
    info=TaoAppSetDestroyRoutine(taoapp,TaoAppDestroyColoringXXX, (void*)mctx); CHKERRQ(info);
  }
  /*
  if (coloring){
  info=PetscObjectReference((PetscObject)coloring); CHKERRQ(info);
  }
  if (mctx->coloring){
       info=ISColoringDestroy(mctx->coloring); CHKERRQ(info);
  }
  */
  mctx->coloring=coloring;  
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetColoring"
/*@C
   TaoAppGetColoring - Set the matrix coloring to be used when computing the
   Hessian by finite differences.

   Collective on TAO_APPLICATION

   Input Parameters:
+  app - the TAO_APPLICATION context
-  coloring - the coloring to be used

   Level: advanced

   Note:
   The gradient evaluation must be set using the routine TaoAppSetGradientRoutine().

 .keywords: TAO_APPLICATION, finite differences, Hessian, coloring, sparse

.seealso: TaoAppSetHessianRoutine(), TaoAppDefaultComputeHessianColor(),
          TaoAppSetGradientRoutine()

@*/
int TaoAppGetColoring(TAO_APPLICATION taoapp, ISColoring *coloring){
  int info;
  TaoAppColorit *mctx;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (coloring){
    info = TaoAppQueryForObject(taoapp,TaoAppColoringXXX,(void**)&mctx); CHKERRQ(info);
    if (mctx){
      *coloring=mctx->coloring;
    } else {
      *coloring=0;
    }
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppAddFiniteDifferences"
int TaoAppAddFiniteDifferences(TAO_APPLICATION taoapp){
  int info;
  PetscFunctionBegin;
  info = TaoAppSetOptionsRoutine(taoapp,TaoAppSetFiniteDifferencesOptions); CHKERRQ(info);
  PetscFunctionReturn(0);
}

int MatTAOMFSetGradient(Mat mat,Vec v,int (*func)(TAO_APPLICATION,Vec,Vec,void *),void *funcctx){
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

