#include "taodaapplication.h"     /*I "taodaapplication.h" I*/
#include "taoapp.h"
#include "petsc.h"
#include "daapp_impl.h"

int DAAPP_COOKIE=0;
static char PPetscDAAppXXX[] = "PPetscDAApp";

int TaoAppDAApp(TAO_APPLICATION, DA_APPLICATION *);
int DAAppDestroy(DA_APPLICATION);
int DAAppExtend(TAO_APPLICATION);

#undef __FUNCT__  
#define __FUNCT__ "DAApplicationCreate"
/* @C
  DAApplicationCreate - Creates a TaoApplication defined
a PETSc Distributed Array (DA) object.  The routines used to evaluate
the objective function, constraint, and derivative information are
designed specifically for these problems.

   Input Parameters:
+  comm - an MPI communicator
.  tda - an array of Distributed Array objects
-  nnda - the number of Distibuted Array objects

   Output Parameters:
.  newdaapplication - the TaoApplication structure

.seealso TaoDAAppSolve(), TaoApplicationCreate(), TaoAppDestroy();

   Level: beginner

.keywords: Application, DA
@ */
int DAApplicationCreate(MPI_Comm comm, DA* tda, PetscInt nnda, TAO_APPLICATION* newdaapplication){
  int info;
  TAO_APPLICATION ttapp;
  PetscFunctionBegin;
  info = TaoApplicationCreate(comm,&ttapp); CHKERRQ(info);
  info = TaoAppSetDAApp(ttapp,tda,nnda); CHKERRQ(info);
  info = PetscInfo(tda[0],"Creating a DA_APPLICATION object.\n"); CHKERRQ(info);
  *newdaapplication=ttapp;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAAppExtend"
int DAAppExtend(TAO_APPLICATION daapplication){
  int info;
  MPI_Comm comm;
  DA_APPLICATION daapp;

  PetscFunctionBegin;
  if (DAAPP_COOKIE==0){
      info=PetscCookieRegister("TAO DA Application",&DAAPP_COOKIE);CHKERRQ(info);
  }
  info=PetscObjectGetComm((PetscObject)daapplication,&comm); CHKERRQ(info);
  info = PetscHeaderCreate(daapp,_p_DA_APPLICATION,PetscInt,DAAPP_COOKIE,-1,"DA Application",comm,DAAppDestroy,0); CHKERRQ(info);
  
  daapp->nda=0;
  daapp->ndamax=PETSCDAAPPMAXGRIDS;
  daapp->currentlevel=0;
  daapp->IsComplementarity=PETSC_FALSE;
  daapp->kspflag=SAME_NONZERO_PATTERN;
  daapp->computedafunction=0;
  daapp->computedagradient=0;
  daapp->computedafunctiongradient=0;
  daapp->computedahessian=0;
  daapp->usrdafctx=0; daapp->usrdafgctx=0; daapp->usrdagctx=0; daapp->usrdahctx=0;
  daapp->computedabounds=0;
  daapp->bounddactx=0;
  daapp->nbeforemonitors=0;
  daapp->naftermonitors=0;
  PetscInt ii;
  info=TaoAppAddObject(daapplication,PPetscDAAppXXX,(void*)daapp,&ii); CHKERRQ(info);
  info=TaoAppSetDestroyRoutine(daapplication,(int(*)(void*))(TaoAppDestroyDAApp), (void*)daapplication); CHKERRQ(info);
  info=TaoAppSetOptionsRoutine(daapplication,DAAppSetOptions); CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetDAApp"
/*@
  TaoAppSetDAApp - Extend the functionality of a TaoApplication by 
adding support for applications defined
a PETSc Distributed Array (DA) object.  The routines used to evaluate
the objective function, constraint, and derivative information are
designed specifically for these problems.

   Input Parameters:
+  daappliation - an existing application object
.  tda - an array of Distributed Array objects
-  nnda - the number of Distibuted Array objects

.seealso TaoDAAppSolve(), TaoApplicationCreate();

   Level: beginner

.keywords: Application, DA
@*/
int TaoAppSetDAApp(TAO_APPLICATION daapplication, DA* tda, PetscInt nnda){
  int info;
  PetscInt i;
  DA_APPLICATION daapp;
  PetscScalar zero=0.0;
  PetscFunctionBegin;
  if (nnda > PETSCDAAPPMAXGRIDS){
    SETERRQ1(1,"Number of grids cannot exceed %d.\n",PETSCDAAPPMAXGRIDS); }
  info = DAAppExtend(daapplication); CHKERRQ(info);
  info = TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  info = PetscInfo(daapplication,"TAO Extending TAO_APPLICATION for DA_APPLICATION object.\n"); CHKERRQ(info);
  info = DAAppSetMatType(daapplication,(const MatType)MATAIJ); CHKERRQ(info);
  daapp->nda=nnda;
  daapp->currentlevel=0;
  for (i=0;i<nnda; i++){
    PetscValidHeaderSpecific(tda[i],DM_COOKIE,2);
    info = PetscObjectReference((PetscObject)(tda[i])); CHKERRQ(info);
    daapp->grid[i].da=tda[i];
    info = DACreateGlobalVector(daapp->grid[i].da,&daapp->grid[i].X); CHKERRQ(info);
    info = VecDuplicate(daapp->grid[i].X,&daapp->grid[i].R); CHKERRQ(info);
    info = VecDuplicate(daapp->grid[i].X,&daapp->grid[i].RHS); CHKERRQ(info);
    daapp->grid[i].XL=0;
    daapp->grid[i].XU=0;
    daapp->grid[i].H=0;
    daapp->grid[i].Interpolate=0;
    daapp->grid[i].CScale=0;
    daapp->grid[i].mgrid=0;
    if (i>0){
      info = DAGetInterpolation(daapp->grid[i-1].da,daapp->grid[i].da,&daapp->grid[i].Interpolate,&daapp->grid[i].CScale); CHKERRQ(info);
    }
  }
  info = VecSet(daapp->grid[0].X, zero); CHKERRQ(info);
  info = TaoAppSetInitialSolutionVec(daapplication,daapp->grid[0].X); CHKERRQ(info);
  info = PetscInfo(daapp,"Create work vectors for DA_APPLICATION object.\n"); CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppDAApp"
int TaoAppDAApp(TAO_APPLICATION daapplication, DA_APPLICATION *daapp){
  int info;
  DA_APPLICATION daapp2;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(daapplication,TAO_APP_COOKIE,1);
  info=TaoAppQueryForObject(daapplication,PPetscDAAppXXX,(void**)&daapp2); CHKERRQ(info);
  if (daapp2){
    PetscValidHeaderSpecific(daapp2,DAAPP_COOKIE,0);
  } else {
    SETERRQ(1,"TAO ERROR: Must call TaoAppSetDAApp() first");
  }
  *daapp=daapp2;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppDestroyDAApp"
int TaoAppDestroyDAApp(TAO_APPLICATION daapplication){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppQueryForObject(daapplication,PPetscDAAppXXX,(void**)&daapp); CHKERRQ(info);
  if (daapp){
    info=DAAppDestroy(daapp); CHKERRQ(info);
    //    info=TaoAppQueryRemoveObject(daapplication,PPetscDAAppXXX); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAAppDestroy"
int DAAppDestroy(DA_APPLICATION daapp){
  PetscInt i,nnda=daapp->nda;
  int info;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(daapp,DAAPP_COOKIE,1);

  if (--((PetscObject)daapp)->refct > 0) PetscFunctionReturn(0);

  info = PetscInfo(daapp,"Destroying work vectors for DA_APPLICATION object.\n"); CHKERRQ(info);
  for (i=0;i<nnda; i++){
    if (daapp->grid[i].coloring){
      info = ISColoringDestroy(daapp->grid[i].coloring); CHKERRQ(info); daapp->grid[i].coloring=0;
    }
    if (daapp->grid[i].H){
      info = MatDestroy(daapp->grid[i].H); CHKERRQ(info);daapp->grid[i].H=0;
    }
    if (daapp->grid[i].XL && daapp->grid[i].XU){
      info = VecDestroy(daapp->grid[i].XL); CHKERRQ(info);daapp->grid[i].XL=0;
      info = VecDestroy(daapp->grid[i].XU); CHKERRQ(info);daapp->grid[i].XU=0;
    }
    if (i>0){
      info = MatDestroy(daapp->grid[i].Interpolate); CHKERRQ(info);daapp->grid[i].Interpolate=0;
      info = VecDestroy(daapp->grid[i].CScale); CHKERRQ(info);daapp->grid[i].CScale=0;
    }
    info = VecDestroy(daapp->grid[i].RHS); CHKERRQ(info);daapp->grid[i].RHS=0;
    info = VecDestroy(daapp->grid[i].R); CHKERRQ(info);daapp->grid[i].R=0;
    info = VecDestroy(daapp->grid[i].X); CHKERRQ(info);daapp->grid[i].X=0;
    info = DADestroy(daapp->grid[i].da); CHKERRQ(info);daapp->grid[i].da=0;
  }
  daapp->nda=0;
  daapp->currentlevel=0;
  PetscHeaderDestroy(daapp); 
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAAppGetDA"
/*@
  DAAppGetDA - Get the DA on a the specified level.
  

   Input Parameters:
+  daapplication - the TAO DAApplication structure
-  n - the level

   Output Parameter:
.  da - address of the pointer to the DA.


.seealso TaoAppSetDAApp();

   Level: intermediate

.keywords: Application, DA
@*/
int DAAppGetDA(TAO_APPLICATION daapplication, PetscInt n, DA *da){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  *da=daapp->grid[n].da;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAAppGetNumberOfDAGrids"
/*@
  DAAppGetNumberOfDAGrids - Get the number of grids specified on the application
  

   Input Parameters:
.  daapplication - the TAO DAApplication structure

   Output Parameter:
.  n - number of DA grids specified

.seealso TaoAppSetDAApp();

   Level: intermediate

.keywords: Application, DA
@*/
int DAAppGetNumberOfDAGrids(TAO_APPLICATION daapplication, PetscInt *n){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  *n=daapp->nda;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAAppGetCurrentLevel"
/*@
  DAAppGetCurrentLevel - Get the number of the current DA grid

   Input Parameters:
.  daapplication - the TAO DAApplication structure

   Output Parameter:
.  n - number of DA grids specified

.seealso TaoAppSetDAApp();

   Level: intermediate

.keywords: Application, DA
@*/
int DAAppGetCurrentLevel(TAO_APPLICATION daapplication, PetscInt *n){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  *n=daapp->currentlevel;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DAAppPetscHessian"
static int DAAppPetscHessian(TAO_APPLICATION daapplication , Vec X , Mat *M, Mat *MPre, MatStructure *flag, void*ctx){
  int     info;
  PetscInt clevel;
  Mat     H=*M;
  DA_APPLICATION daapp=(DA_APPLICATION)ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(daapplication,TAO_APP_COOKIE,1);
  PetscValidHeaderSpecific(daapp,DAAPP_COOKIE,1);
  PetscStackPush("TAO User DA Hessian");
  clevel=daapp->currentlevel;
  
  //    daapp->FDGrad=daapp->grid[i].RHS;   /* Needed for finite differences gradient vector */
  info = PetscInfo1(daapp,"TAO hessian evaluation at DA_APPLICATION object, level %d.\n",clevel); CHKERRQ(info);
  info = (*daapp->computedahessian)(daapplication,daapp->grid[clevel].da,X,H,daapp->usrdahctx);CHKERRQ(info);
  *flag=daapp->kspflag;

  PetscStackPop;

  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DAAppSetHessianMat"
/*@
  DAAppSetHessianMat - Creates the matrices used to store the Hessian at finer grid levels

   Collective on TAO_APPLICATION

   Input Parameters:
.  daapplication - the DA Application object

   Level: advanced


.keywords:  DA, hessian

.seealso: DAAppSetHessianRoutine(), DAAppSetMatType();

@*/
int DAAppSetHessianMat(TAO_APPLICATION daapplication){
  int info;
  PetscInt i,nnda;
  DA_APPLICATION daapp;

  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  nnda=daapp->nda;
  for (i=0;i<nnda; i++){
    if (daapp->grid[i].H==0){
      info = DAGetMatrix(daapp->grid[i].da,daapp->HessianMatrixType,&daapp->grid[i].H); CHKERRQ(info);
      info = MatSetOption(daapp->grid[i].H,MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERRQ(info);
      info = MatSetOption(daapp->grid[i].H,MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE); CHKERRQ(info);
      info = MatSetOption(daapp->grid[i].H,MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(info);
      info = PetscInfo1(daapp->grid[i].da,"TAO create Hessian for DA_APPLICATION object, level %d.\n",i); CHKERRQ(info);
      info = DAGetColoring(daapp->grid[i].da,IS_COLORING_GLOBAL,daapp->HessianMatrixType,&(daapp->grid[i].coloring));CHKERRQ(info);
    }
  }
  info = TaoAppSetHessianMat(daapplication,daapp->grid[0].H,daapp->grid[0].H);CHKERRQ(info);
  
  PetscFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "DAAppSetHessianRoutine"
/*@C
  DAAppSetHessianRoutine - Set a routine that will evaluate the hessian
  on the given DA at the given point.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
.  hess - the function pointer for the hessian evaluation routine
-  ctx - the hessian context

   Calling sequence of hess:
$     hess(TAO_APPLICATION daapplication,DA da, Vec x,Mat H,void *ctx);

+  daapplication - the TAO_APPLICATION context
.  da - the Distributed Array
.  x - input vector
.  H - hessian matrix
-  ctx - user-defined hessian context set from DAAppSetHessianRoutine()

   Level: beginner

   Options Database Key:
.  -tao_view_hessian - view the hessian after each evaluation using PETSC_VIEWER_STDOUT_WORLD

.keywords:  DA, hessian

.seealso: DAAppSetObjectiveAndGradientRoutine();
@*/
int DAAppSetHessianRoutine(TAO_APPLICATION daapplication, int (*hess)(TAO_APPLICATION,DA,Vec,Mat,void*),void *ctx){
  int info;
  DA_APPLICATION daapp;

  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  daapp->computedahessian=hess;
  daapp->usrdahctx=ctx;
  if (hess){
    info = TaoAppSetHessianRoutine(daapplication,DAAppPetscHessian,(void*)daapp);CHKERRQ(info);
  } else {
    info = TaoAppSetHessianRoutine(daapplication,0,0);CHKERRQ(info);
  }
  info = DAAppSetHessianMat(daapplication);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAAppPetscFunctionGradient"
static int DAAppPetscFunctionGradient(TAO_APPLICATION daapplication , Vec X ,double *ff, Vec G, void*ctx){
  int info;
  DA_APPLICATION daapp=(DA_APPLICATION)ctx;
  PetscFunctionBegin;

  PetscValidHeaderSpecific(daapplication,TAO_APP_COOKIE,1);
  PetscValidHeaderSpecific(daapp,DAAPP_COOKIE,1);
  PetscStackPush("TAO User DA Objective and gradient function");
  info = (*daapp->computedafunctiongradient)(daapplication,daapp->grid[daapp->currentlevel].da,X,ff,G,daapp->usrdafgctx);
  CHKERRQ(info);
  info = PetscInfo1(daapplication,"TAO function evaluation at DA_APPLICATION object, level %d.\n",daapp->currentlevel); CHKERRQ(info);
  PetscStackPop;

  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DAAppSetObjectiveAndGradientRoutine"
/*@C
  DAAppSetObjectiveAndGradientRoutine - Set a routine that will evaluate the objective and
  gradient functions on the given DA at the given point.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
.  grad - the function pointer for the gradient evaluation routine
-  ctx - the function-gradient context

   Calling sequence of funcgrad:
$     funcgrad(TAO_APPLICATION daapplication,DA da, Vec x,double *f,Vec g,void *ctx);

+  daapplication - the TAO_APPLICATION daapplication context
.  da - the Distributed Array
.  x - input vector
.  f - objective value 
.  g - gradient vector 
-  ctx - user defined function-gradient context set from DAAppSetObjectiveAndGradientRoutine()

   Fortran Note:
   If your Fortran compiler does not recognize symbols over 31 characters in length, then
   use the identical routine with the shortened name DAAppSetObjectiveAndGradientRou()

   Level: beginner

   Options Database Key:
.  -tao_view_gradient - view the gradient after each evaluation using PETSC_VIEWER_STDOUT_SELF

.keywords: DA, Gradient, Objective Function

.seealso: DAAppSetObjectiveRoutine(), DAAppSetGradientRoutine();

@*/
int DAAppSetObjectiveAndGradientRoutine(TAO_APPLICATION daapplication, int (*funcgrad)(TAO_APPLICATION,DA,Vec,double*,Vec, void*),void *ctx){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (funcgrad){ 
    info=TaoAppSetObjectiveAndGradientRoutine(daapplication,
					      DAAppPetscFunctionGradient,(void*)daapp);CHKERRQ(info);
  } else {
    info=TaoAppSetObjectiveAndGradientRoutine(daapplication,0,0);CHKERRQ(info);
  }
  info=TaoAppSetInitialSolutionVec(daapplication,daapp->grid[0].X);CHKERRQ(info);
  daapp->computedafunctiongradient=funcgrad;
  daapp->usrdafgctx=ctx;
  info = PetscInfo(daapp,"Set objective function pointer for DA_APPLICATION object.\n"); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAAppPetscGradient"
static int DAAppPetscGradient(TAO_APPLICATION daapplication , Vec X ,Vec G, void*ctx){
  int info;
  DA_APPLICATION daapp=(DA_APPLICATION)ctx;
  PetscFunctionBegin;

  PetscValidHeaderSpecific(daapplication,TAO_APP_COOKIE,1);
  PetscValidHeaderSpecific(daapp,DAAPP_COOKIE,1);
  PetscStackPush("TAO User DA Gradient Evaluation");
  info = (*daapp->computedagradient)(daapplication,daapp->grid[daapp->currentlevel].da,X,G,daapp->usrdafgctx);
  CHKERRQ(info);
  info = PetscInfo1(daapplication,"TAO gradient evaluation at DA_APPLICATION object, level %d.\n",daapp->currentlevel); CHKERRQ(info);
  PetscStackPop;

  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DAAppSetGradientRoutine"
/*@C
  DAAppSetGradientRoutine - Set a routine that will evaluate the gradient 
  function on the given DA at the given point.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
.  grad - the function pointer for the gradient evaluation routine
-  ctx - the gradient context

   Calling sequence of grad:
$     grad(TAO_APPLICATION daapplication,DA da, Vec x,Vec g,void *ctx);

+  daapplication - the TAO_APPLICATION context
.  da - the Distributed Array
.  x - input vector
.  g - gradient vector 
-  ctx - user defined gradient context set from DAAppSetGradientRoutine()

   Level: intermediate

   Options Database Key:
.  -tao_view_gradient - view the gradient after each evaluation using PETSC_VIEWER_STDOUT_SELF

.keywords:  DA, gradient

.seealso: DAAppSetObjectiveRoutine(), DAAppSetObjectiveAndGradientRoutine();

@*/
int DAAppSetGradientRoutine(TAO_APPLICATION daapplication, int (*grad)(TAO_APPLICATION,DA,Vec,Vec, void*),void *ctx){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (grad){ 
    info=TaoAppSetGradientRoutine(daapplication, DAAppPetscGradient, (void*)daapp);CHKERRQ(info);
  } else {
    info=TaoAppSetGradientRoutine(daapplication,0,0);CHKERRQ(info);
  }
  daapp->computedagradient=grad;
  daapp->usrdagctx=ctx;
  info = PetscInfo(daapp,"Set gradient function pointer for DA_APPLICATION object.\n"); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAAppPetscObjectiveFunction"
static int DAAppPetscObjectiveFunction(TAO_APPLICATION daapplication , Vec X ,double* f, void*ctx){
  int info;
  DA_APPLICATION daapp=(DA_APPLICATION)ctx;
  PetscFunctionBegin;

  PetscValidHeaderSpecific(daapplication,TAO_APP_COOKIE,1);
  PetscValidHeaderSpecific(daapp,DAAPP_COOKIE,1);
  PetscStackPush("TAO User DA Objective function");
  info = (*daapp->computedafunction)(daapplication,daapp->grid[daapp->currentlevel].da,X,f,daapp->usrdafgctx);
  CHKERRQ(info);
  info = PetscInfo1(daapp,"TAO gradient evaluation at DA_APPLICATION object, level %d.\n",daapp->currentlevel); CHKERRQ(info);
  PetscStackPop;

  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DAAppSetObjectiveRoutine"
/*@C
  DAAppSetObjectiveRoutine - Set a routine that will evaluate the objective 
  function on the given DA at the given point.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
.  func - the function pointer for the objecive evaluation routine
-  ctx - the monitor context

   Calling sequence of func:
$     func(TAO_APPLICATION daapplication,DA da,Vec x,double *f,void *ctx);

+  daapplication - the TAO_APPLICATION context
.  da - the Distributed Array
.  x - input vector
.  f - application sets equal to the function value 
-  ctx - user-defined function context set from DAAppSetObjectiveRoutine()

   Level: beginner

.keywords:  DA, objective

.seealso: DAAppSetObjectiveAndGradientRoutine();

@*/
int DAAppSetObjectiveRoutine(TAO_APPLICATION daapplication, int (*func)(TAO_APPLICATION,DA,Vec,double*, void*),void *ctx){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (func){
    info=TaoAppSetObjectiveRoutine(daapplication, DAAppPetscObjectiveFunction, (void*)daapp);CHKERRQ(info);
  } else {
    info=TaoAppSetObjectiveRoutine(daapplication,0,0);CHKERRQ(info);
  }
  info=TaoAppSetInitialSolutionVec(daapplication, daapp->grid[0].X);CHKERRQ(info);
  daapp->computedafunction=func;
  daapp->usrdafctx=ctx;
  info = PetscInfo(daapp,"Set objective function pointer for DA_APPLICATION object.\n"); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAApplicationEvaluateVariableBounds"
int DAApplicationEvaluateVariableBounds(TAO_APPLICATION daapplication, Vec XL, Vec XU, void*ctx){
  DA_APPLICATION daapp=(DA_APPLICATION)ctx;
  PetscInt level=daapp->currentlevel;
  int  info;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(daapplication,TAO_APP_COOKIE,1);
  info = PetscInfo1(daapp->grid[level].da,"Compute variable bounds in level %d of DA_APPLICATION object.\n",level); CHKERRQ(info);
  info = PetscObjectReference((PetscObject)XL);
  info = PetscObjectReference((PetscObject)XU);
  daapp->grid[level].XL=XL;
  daapp->grid[level].XU=XU;
  if (daapp->computedabounds){
    info = (*daapp->computedabounds)(daapplication,daapp->grid[level].da,XL,XU,daapp->bounddactx); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DAAppSetVariableBoundsRoutine"
/*@C
  DAAppSetVariableBoundsRoutine - Set a routine that will evaluate the bounds
  on the variables.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
.  bound - the function pointer for the bound evaluation routine
-  ctx - the monitor context

   Calling sequence of bound:
$     bound(TAO_APPLICATION daapplication, DA da,Vec xl, Vec xu, void *ctx);

+  daapplication - the TAO_APPLICATION context
.  da - the Distributed Array
.  xl - vector of upper bounds
.  xu - vector of lower bounds
-  ctx - user-defined monitor context set from DAAppSetVariableBoundsRoutine()

   Level: beginner

.keywords:  DA, bounds

.seealso: TaoDAAppSolve();

@*/
int DAAppSetVariableBoundsRoutine(TAO_APPLICATION daapplication, int (*bound)(TAO_APPLICATION,DA,Vec,Vec, void*),void *ctx){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info = TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  info=TaoAppSetVariableBoundsRoutine(daapplication,DAApplicationEvaluateVariableBounds,(void*)daapp);CHKERRQ(info);
  daapp->computedabounds=bound;
  daapp->bounddactx=ctx;
  info = PetscInfo(daapp,"Set variable bounds in DA_APPLICATION object.\n"); CHKERRQ(info);
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "DAAppSetConstraintRoutine"
/*@C
  DAAppSetConstraintRoutine - Set a routine that will evaluate the constraint
  function on the given DA at the given point.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the TAO Application object
.  grad - the function pointer for the gradient evaluation routine
-  ctx - the gradient context

   Calling sequence of grad:
$     f(TAO_APPLICATION daapplication,DA da, Vec x,Vec r,void *ctx);

+  daapplication - the DA_APPLICATION context
.  da - the Distributed Array
.  x - input vector
.  r - constraint vector 
-  ctx - user defined gradient context set from DAAppSetGradientRoutine()

   Level: intermediate

   Options Database Key:
.  -tao_view - view the gradient after each evaluation using PETSC_VIEWER_STDOUT_SELF

.keywords:  DA, gradient

.seealso: DAAppSetJacobianRoutine();
@*/
int DAAppSetConstraintRoutine(TAO_APPLICATION daapplication, int (*f)(TAO_APPLICATION,DA,Vec,Vec, void*),void *ctx){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (f){ 
    info=TaoAppSetConstraintRoutine(daapplication, DAAppPetscGradient, (void*)daapp);CHKERRQ(info);
  } else {
    info=TaoAppSetConstraintRoutine(daapplication,0,0);CHKERRQ(info);
  }
  daapp->IsComplementarity=PETSC_TRUE;
  info=DAAppSetGradientRoutine(daapplication,f,ctx); CHKERRQ(info);
  daapp->computedagradient=f;
  daapp->usrdagctx=ctx;
  daapp->usrdafgctx=ctx;

  info = PetscInfo(daapp,"Set constraint function pointer for DA_APPLICATION object.\n"); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DAAppSetJacobianRoutine"
/*@C
  DAAppSetJacobianRoutine - Set a routine that will evaluate the Jacobian
  on the given DA at the given point.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
.  jac - the function pointer for the Jacobian evaluation routine
-  ctx - the jacobian context

   Calling sequence of hess:
$     jac(TAO_APPLICATION daapplication, DA da, Vec x,Mat J,void *ctx);

+  daapplication - the TAO_APPLICATION context
.  da - the Distributed Array
.  x - input vector
.  J - Jacobian matrix
-  ctx - user-defined hessian context set from DAAppSetHessianRoutine()

   Level: intermediate

   Options Database Key:
.  -tao_view_jacobian - view the jacobian after each evaluation using PETSC_VIEWER_STDOUT_WORLD

.keywords:  DA, hessian

.seealso: DAAppSetConstraintRoutine();

@*/
int DAAppSetJacobianRoutine(TAO_APPLICATION daapplication, int (*jac)(TAO_APPLICATION,DA,Vec,Mat,void*),void *ctx){
  int info;
  PetscInt i,nnda;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info = TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (!jac){ SETERRQ(1,"No DA Hessian routine provided\n");}
  daapp->computedahessian=jac;
  daapp->usrdahctx=ctx;
  nnda=daapp->nda;;
  for (i=0;i<nnda; i++){
    if (daapp->grid[i].H==0){
      info = DAGetMatrix(daapp->grid[i].da,daapp->HessianMatrixType,&daapp->grid[i].H); CHKERRQ(info);
      info = MatSetOption(daapp->grid[i].H,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(info);
      info = MatSetOption(daapp->grid[i].H,MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE); CHKERRQ(info);
      info = PetscInfo1(daapp->grid[i].da,"TAO create Jacobian for DA_APPLICATION object, level %d.\n",i); CHKERRQ(info);
    }
  }
  
  info = TaoAppSetJacobianRoutine(daapplication,DAAppPetscHessian,(void*)daapp);CHKERRQ(info);
  info = TaoAppSetJacobianMat(daapplication,daapp->grid[0].H,daapp->grid[0].H);CHKERRQ(info);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppSetBeforeMonitor"
/*@C
  DAAppSetBeforeMonitor - Set a routine that will be called before the
  optimization on each grid.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
.  beforemonitor - a monitor routine called before solving the application on each DA 
-  ctx - the monitor context

   Calling sequence of monitor:
$     beforemonitor(TAO_APPLICATION daapplication, DA da,int level, void *ctx);

+  daapplication - this TAO_APPLICATION context
.  da - the Distributed Array
.  level - the grid level that will be solved next (level 0 is coarsest)
-  ctx - user-defined function context set from DAAppSetBeforeMonitor()

Note:
   These monitors are different from the monitors that can be called after
   each iteration of the optimization algorithm.

Note:
   The beforemonitor and aftermonitor are good for setting up and destroying the application data.

   Level: intermediate

.keywords:  DA, monitor

.seealso: DAAppSetAfterMonitor(), TaoSetMonitor(), TaoAppSetMonitor();

@*/
int DAAppSetBeforeMonitor(TAO_APPLICATION daapplication, int (*beforemonitor)(TAO_APPLICATION,DA,PetscInt, void*), void *ctx){
  PetscInt n;
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (beforemonitor){
    n=daapp->nbeforemonitors;
    if (n>=MAX_DAAP_MONITORS){
      SETERRQ(1,"TAO ERROR: Too many beforemonitors set in DAAPP");
    }
    daapp->beforemonitor[n]=beforemonitor;
    daapp->beforemonitorctx[n]=ctx;
    daapp->nbeforemonitors++;
    info = PetscInfo(daapplication,"TAO: Set user beforemonitor for DA_APPLICATION object.\n"); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppSetAfterMonitor"
/*@C
  DAAppSetAfterMonitor - Set a routine that will be called after the
  optimization on each grid.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
.  aftermonitor - a monitor routine called after solving the application on each DA 
-  ctx - the monitor context

   Calling sequence of monitor:
$     aftermonitor(TAO_APPLICATION daapplication, DA da, int level, void *ctx);

+  daapplication - this TAO_APPLICATION context
.  da - the Distributed Array
.  level - the grid level that will be solved next (level 0 is coarsest)
-  ctx - user-defined function context set from DAAppSetAfterMonitor()

Note:
   These monitors are different from the monitors that can be called after
   each iteration of the optimization algorithm.

Note:
   The beforemonitor and aftermonitor are good for setting up and destroying the application data.

   Level: intermediate

.keywords:  DA, monitor

.seealso: DAAppSetBeforeMonitor(), TaoSetMonitor(), TaoAppSetMonitor();

@*/
int DAAppSetAfterMonitor(TAO_APPLICATION daapplication, int (*aftermonitor)(TAO_APPLICATION,DA,PetscInt, void*), void *ctx){
  int info;
  PetscInt n;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (aftermonitor){
    n=daapp->naftermonitors;
    if (n>=MAX_DAAP_MONITORS){
      SETERRQ(1,"TAO ERROR: Too many aftermonitors set in DAAPP");
    }
    daapp->aftermonitor[n]=aftermonitor;
    daapp->aftermonitorctx[n]=ctx;
    daapp->naftermonitors++;
    info = PetscInfo(daapp,"TAO: Set user aftermonitor for DA_APPLICATION object.\n"); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DAAppGetSolution"
/*@
  DAAppGetSolution - Sets a pointer to an existing vector that contains
  the solution.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
-  level - the DA on which the current solution is desired (Level 0 is coarsest)

   Output Parameters:
.  X - the current solution

   Level: beginner

.seealso: DAAppSetBeforeMonitor(), DAAppSetAfterMonitor();

.keywords: Solution, DA
@*/
int DAAppGetSolution(TAO_APPLICATION daapplication, PetscInt level, Vec *X){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (level<0 || level>=daapp->nda){
    SETERRQ(1,"No Solution available.  This grid does not exist");
  } else if (X && daapp && daapp->grid && level<daapp->nda){
    *X=daapp->grid[level].X;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DAAppGetHessianMat"
/*@C
  DAAppGetHessianMat - Sets a pointer to an existing matrix that contains
  the Hessian at the specified level.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
-  level - the DA on which the current Hessian is desired (Level 0 is coarsest)

   Output Parameters:
.  H - the current Hessian

   Level: intermediate

.seealso: DAAppSetHessianRoutine;

.keywords: Hessian, DA
@*/
int DAAppGetHessianMat(TAO_APPLICATION daapplication, PetscInt level, Mat *H){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (level<0 || level>=daapp->nda){
    SETERRQ(1,"No Hessian available.  This grid does not exist");
  } else if (H && daapp && daapp->grid && level<daapp->nda){
    *H=daapp->grid[level].H;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppGetInterpolationMatrix"
/*@
  DAAppGetInterpolationMatrix - Sets pointers to existing vectors
  containg upper and lower bounds on the variables.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
-  level - the DA on which the current solution is desired (Level 0 is coarsest)

   Output Parameters:
+  Interpolate - the interpolating matrix
-  CScale - the scaling matrix for restriction Matrix

   Level: intermediate

.seealso: DAAppSetBeforeMonitor(), DAAppSetAfterMonitor();

.keywords: Interpolate, DA
@*/
int DAAppGetInterpolationMatrix(TAO_APPLICATION daapplication, PetscInt level, Mat *Interpolate, Vec *CScale){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if ( !daapp || level<1 || level>=daapp->nda){
    SETERRQ1(1,"No Interpolator available for level %d.  This grid does not exist",level);
  } 
  if (Interpolate){
    *Interpolate=daapp->grid[level].Interpolate;
  }
  if (CScale){
    *CScale=daapp->grid[level].CScale;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DAAppGetVariableBounds"
/*@
  DAAppGetVariableBounds - Sets pointers to existing vectors
  containg upper and lower bounds on the variables.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
-  level - the DA on which the current solution is desired (Level 0 is coarsest)

   Output Parameters:
+  XL - the current solution
-  XU - the current solution

   Level: intermediate

.seealso: DAAppSetBeforeMonitor(), DAAppSetAfterMonitor();

.keywords: Solution, DA
@*/
int DAAppGetVariableBounds(TAO_APPLICATION daapplication, PetscInt level, Vec *XL, Vec *XU){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if ( !daapp || level<0 || level>=daapp->nda){
    SETERRQ(1,"No Solution available.  This grid does not exist.");
  } 
  if (XL){
    *XL=daapp->grid[level].XL;
  }
  if (XU){
    *XU=daapp->grid[level].XU;
  }
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DAAppSetInitialSolution"
/*@
  DAAppSetInitialSolution - Sets the initial solution for
  the grid application.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the DA Application object
-  X0 - the initial solution vector

   Level: beginner

.seealso: DAAppGetSolution();

.keywords: Initial Solution, DA
@*/
int DAAppSetInitialSolution(TAO_APPLICATION daapplication, Vec X0){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(X0,VEC_COOKIE,2);
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  info=VecCopy(X0,daapp->grid[0].X); CHKERRQ(info);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "DAAppSetMatStructure"
/*
  DAAppSetMatStructure - Set the MatStructure flag to be used in KSPSetOperators:

  Collective on TAO_APPLICATION
  
  Input Parameters:
+  daapplication - the DA Application object
-  flag - indicates the KSP flag

   Level: intermediate

   Note:
   The default value is SAME_NONZERO_PATTERN

.seealso: KSPSetOperators(), TaoAppGetKSP();

.keywords: Linear Solver, Multigrid, DA, KSP

@*/
int DAAppSetMatStructure(TAO_APPLICATION daapplication, MatStructure pflag){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;

  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  daapp->kspflag=pflag;
  info = PetscInfo(daapp,"Set MatStructure flag in TAO solver.\n"); CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppSetMatType"
/*@
  DAAppSetMatType - set the matrix type to be used for the Hessian matrix

  Collective on TAO_APPLICATION
  
  Input Parameters:
+  daapplication - the DA Application object
-  mattype - the matrix type

   Level: intermediate

   Options Database Key:
.  -tao_da_mattype - select matrix type

   Note:
   The default value is MATMPIAIJ.

   Note:
   Call before DAAppSetHessianRoutine()

.seealso: DAAppSetMatStructure(),  DAAppSetHessianRoutine(), DAGetMatrix();

.keywords: DA, Hessian

@*/
int DAAppSetMatType(TAO_APPLICATION daapplication, const MatType mattype){
  int info;
  DA_APPLICATION daapp;
  PetscFunctionBegin;

  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  info=PetscStrcpy(daapp->HessianMatrixType,mattype); CHKERRQ(info);

  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DAAppSetOptions"
/*@
  DAAppSetOptions - Sets various parameters to be used in this application
  and the TAO solver.

   Collective on TAO_APPLICATION

   Input Parameters:
.  daapplication - the DA Application object

   Level: intermediate

Note:  This routine will be called by the TaoAppSetFromOptions().

.keywords:  DA, options

.seealso: TaoDAAppSolve();

@*/
int DAAppSetOptions(TAO_APPLICATION daapplication){
  int info;
  //  char *prefix;
  PetscTruth flg1=PETSC_FALSE,flg2=PETSC_FALSE;
  DA_APPLICATION daapp;
  //  MPI_Comm comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(daapplication,TAO_APP_COOKIE,1);
  info = PetscInfo(daapplication,"TaoDAAppSetOptions(): Reading command line for options\n"); CHKERRQ(info);
  info = TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  /*
  info = PetscObjectGetOptionsPrefix((PetscObject)daapplication,&prefix); CHKERRQ(info);
  info = PetscObjectGetComm((PetscObject)daapplication,&comm); CHKERRQ(info);
  info = PetscOptionsBegin(comm,prefix,"TAO PETSC DA APPLICATIONS","solver");CHKERRQ(info);
  */
  info = DAAppSetMultiGridOptions(daapplication);CHKERRQ(info);


  flg1=PETSC_FALSE,flg2=PETSC_FALSE;
  info = PetscOptionsTruth("-tao_da_samepreconditioner","Use same preconditioner in KSP","DAAppSetMatStructure",
			     PETSC_FALSE,&flg2,&flg1);CHKERRQ(info);
  if (flg1 && flg2) {
    info = DAAppSetMatStructure(daapplication, SAME_PRECONDITIONER);CHKERRQ(info); 
  } else if (flg1){
    info = DAAppSetMatStructure(daapplication, SAME_NONZERO_PATTERN);CHKERRQ(info); 
  }

  info = PetscOptionsString("-tao_da_mattype","the hessian should have matrix type","DAAppSetMatType",
			    daapp->HessianMatrixType,daapp->HessianMatrixType,18,&flg1);CHKERRQ(info);

  /*
  info = PetscOptionsEnd();CHKERRQ(info);
  */
  PetscFunctionReturn(0);
}

