#include "taodaapplication.h"     /* taodaapplication.h */
#include "taoapp.h"

typedef struct {
  /* Function Gradient Evaluation over single element  */
  int  (*computeelementfunctiongradient)(PetscInt[2], PetscScalar[4],double*,PetscScalar[4],void*); 
  void *elementfgctx;
  int elementfgflops;
} TaoDA2D1DOFFGCtx;


#undef __FUNCT__
#define __FUNCT__ "TaoDA2dLoopFunctionGradient"
static int TaoDA2dLoopFunctionGradient(TAO_APPLICATION daapp, DA da, Vec X, double *f, Vec G, void * ctx) {

  TaoDA2D1DOFFGCtx* fgctx = (TaoDA2D1DOFFGCtx*) ctx;
  Vec localX, localG;
  MPI_Comm  comm;
  int info;
  PetscInt i, j;
  PetscInt coor[2];
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  double floc = 0.0, smallf;
  PetscScalar **x, **g;
  PetscScalar zero = 0.0;
  PetscScalar smallX[4], smallG[4];

  PetscFunctionBegin;
  info = DAGetLocalVector(da, &localX); CHKERRQ(info);
  info = DAGetLocalVector(da, &localG); CHKERRQ(info);
  info = VecSet(G, zero); CHKERRQ(info);
  info = VecSet(localG, zero); CHKERRQ(info);

  info = DAGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(info);

  info = DAVecGetArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecGetArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DAGetCorners(da, &xs, &ys, PETSC_NULL, &xm, &ym, PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(da, &gxs, &gys, PETSC_NULL, &gxm, &gym, PETSC_NULL); CHKERRQ(info);

  xe = gxs + gxm - 1;
  ye = gys + gym - 1;
  for (j = ys; j < ye; j++) {
    for (i = xs; i < xe; i++) {

        smallX[0] = x[j][i];
        smallX[1] = x[j][i+1];
        smallX[2] = x[j+1][i];
        smallX[3] = x[j+1][i+1];
        coor[0] = i; coor[1] = j;

        info = fgctx->computeelementfunctiongradient(coor,smallX,&smallf,smallG,fgctx->elementfgctx);

        floc += smallf;

        g[j][i] += smallG[0];
        g[j][i+1] += smallG[1];
        g[j+1][i] += smallG[2];
        g[j+1][i+1] += smallG[3];
    }
  }

  info = PetscLogFlops((ye-ys)*(xe-xs)*(fgctx->elementfgflops + 5)); CHKERRQ(info);

  info = PetscObjectGetComm((PetscObject)X,&comm); CHKERRQ(info);
  info = MPI_Allreduce(&floc, f, 1, MPI_DOUBLE, MPI_SUM, comm); CHKERRQ(info);

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecRestoreArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DALocalToGlobalBegin(da, localG, G); CHKERRQ(info);
  info = DALocalToGlobalEnd(da, localG, G); CHKERRQ(info);

  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localG); CHKERRQ(info);

  PetscFunctionReturn(0);
} /* TaoDA2dLoopFunctionGradient  */

#undef __FUNCT__
#define __FUNCT__ "DAAppSetElementObjectiveAndGradientRoutine"
/*@C
   DAAppSetElementObjectiveAndGradientRoutine - Set routine that evaluates the
   local part of a function on a 2-dimensional DA with 1 degree of freedom. 

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapp - the TAO_APPLICATION solver context
.  funcgrad - local function gradient routine
.  flops - the number of flops done performed in the funcgrad routine
-  fgctx - [optional] user-defined context for private data for the evaluation.

   Calling sequence of funcgrad:
$     int funcgrad(int coordinates[2], PetscScalar x[4], double *f, PetscScalar g[4], void* ctx)

+    coord - the global coordinates [i j] in each direction of the DA
.    x - the variables on the DA ( da[j][i], da[j][j+1], da[j+1][i], da[j+1][i+1] ) (bottom left, bottom right, top left, top right)
.    f - the local function value
.    g - the gradient of this local function for with respect to each variable
-    ctx - user defined context

   Fortran Note:
   If your Fortran compiler does not recognize symbols over 31 characters in length, then
   use the identical routine with the shortened name DAAppSetElementObjectiveAndGrad()
   

   Level: intermediate

.keywords: DA, Object Function, Gradient

.seealso: DAAppSetObjectiveAndGradientRoutine();
@*/
int DAAppSetElementObjectiveAndGradientRoutine(TAO_APPLICATION daapplication, int (*funcgrad)(PetscInt[2],PetscScalar[4],double*,PetscScalar[4],void*),
					       PetscInt flops, void *ctx){
  int info;
  PetscInt n,i;
  PetscInt dim,dof,s;
  DAStencilType st;
  TaoDA2D1DOFFGCtx *fgctx;
  DA da;

  PetscFunctionBegin;
  info=DAAppGetNumberOfDAGrids(daapplication,&n); CHKERRQ(info);
  for (i=0;i<n;i++){
    info = DAAppGetDA(daapplication, i, &da); CHKERRQ(info);
    info = DAGetInfo(da,&dim,0,0,0,0,0,0,&dof,&s,0,&st); CHKERRQ(info);
    if (dim!=2){
      SETERRQ(1,"TAO DA ERROR: DA must have dimension of 2");}
    if (dof!=1){
      SETERRQ(1,"TAO DA ERROR: DA must have exactly 1 degree of freedom per node");}
    if (s!=1){
      SETERRQ(1,"TAO DA ERROR: DA stencil width must equal 1"); }
    if (st!=DA_STENCIL_BOX){
      SETERRQ(1,"TAO DA ERROR: DA stencil must be DA_STENCIL_BOX");}
  }
  PetscNew(TaoDA2D1DOFFGCtx,&fgctx);
  fgctx->computeelementfunctiongradient=funcgrad;
  fgctx->elementfgctx=ctx;
  fgctx->elementfgflops = flops;
  info = DAAppSetObjectiveAndGradientRoutine(daapplication, TaoDA2dLoopFunctionGradient, (void*)fgctx); CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(daapplication,TaoApplicationFreeMemory, (void*)fgctx); CHKERRQ(info);
  info = PetscInfo(daapplication,"Set objective function pointer for TAO_APPLICATION object.\n"); CHKERRQ(info);
  PetscFunctionReturn(0);
}


typedef struct {
  /* Hesian over single element  */
  int  (*computeelementhessian)(PetscInt[2], PetscScalar[4],PetscScalar[4][4],void*); 
  void *elementhctx;
  PetscInt elementhflops;
} TaoDA2D1DOFHCtx;


#undef __FUNCT__
#define __FUNCT__ "TaoDA2dLoopHessian"
static int TaoDA2dLoopHessian(TAO_APPLICATION daapp, DA da, Vec X, Mat A, void* ctx) {

  TaoDA2D1DOFHCtx* hctx = (TaoDA2D1DOFHCtx*)ctx;
  Vec localX;
  int info;
  PetscInt i,j,coor[2];
  PetscInt xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye,ind[4];
  PetscScalar smallX[4];
  PetscScalar smallH[4][4];
  PetscScalar **x;
  PetscTruth assembled;

  PetscFunctionBegin;

  info = DAGetLocalVector(da, &localX); CHKERRQ(info);
  info = MatAssembled(A,&assembled); CHKERRQ(info);
  if (assembled) {
    info = MatZeroEntries(A); CHKERRQ(info);
  }

  info = DAGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(info);

  info = DAVecGetArray(da, localX, (void**)&x); CHKERRQ(info);

  info = DAGetCorners(da, &xs, &ys, PETSC_NULL, &xm, &ym, PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(da, &gxs, &gys, PETSC_NULL, &gxm, &gym, PETSC_NULL); CHKERRQ(info);

  xe = gxs + gxm - 1;
  ye = gys + gym - 1;
  for (j = ys; j < ye; j++) {
    for (i = xs; i < xe; i++) {

      smallX[0] = x[j][i];
      smallX[1] = x[j][i+1];
      smallX[2] = x[j+1][i];
      smallX[3] = x[j+1][i+1];
      coor[0] = i; coor[1] = j;

      info = hctx->computeelementhessian(coor,smallX,smallH,hctx->elementhctx); CHKERRQ(info);

      ind[0] = (j-gys) * gxm + (i-gxs);
      ind[1] = ind[0] + 1;
      ind[2] = ind[0] + gxm;
      ind[3] = ind[2] + 1;
      info = MatSetValuesLocal(A,4,ind,4,ind,(PetscScalar*)smallH[0],ADD_VALUES); CHKERRQ(info);
    }
  }

  info = PetscLogFlops((ye-ys)*(xe-xs)*hctx->elementhflops); CHKERRQ(info);

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);

  info = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(info);


  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);

  PetscFunctionReturn(0);
} /* Compute Hessian */


#undef __FUNCT__
#define __FUNCT__ "DAAppSetElementHessianRoutine"
/*@C
   DAAppSetElementHessianRoutine - Set routine that evaluates the
   local part of a Hessian matrix on a 2-dimensional DA with 1 degree of freedom. 

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapp - the TAO_APPLICATION solver context
.  hess - local function gradient routine
.  flops - the number of flops done performed in the hess routine
-  fgctx - [optional] user-defined context for private data for the evaluation.

   Calling sequence of funcgrad:
$     int hess(int coordinates[2], PetscScalar x[4], PetscScalar h[4][4], void* ctx)

+    coord - the global coordinates [i j] in each direction of the DA
.    x - the variables on the DA ( da[j][i], da[j][j+1], da[j+1][i], da[j+1][i+1] ) (bottom left, bottom right, top left, top right)
.    h - the hessian of this local function for with respect to each variable
-    ctx - user defined context

   Level: intermediate

.keywords: DA, gradient

.seealso: DAAppSetHessianRoutine()
@*/
int DAAppSetElementHessianRoutine(TAO_APPLICATION daapplication, int (*hess)(PetscInt[2],PetscScalar[4],PetscScalar[4][4],void*), PetscInt flops, void *ctx){
  int info;
  PetscInt n,i;
  PetscInt dim,dof,s;
  DAStencilType st;
  TaoDA2D1DOFHCtx *hctx;
  DA da;

  PetscFunctionBegin;
  info=DAAppGetNumberOfDAGrids(daapplication,&n); CHKERRQ(info);
  for (i=0;i<n;i++){
    info = DAAppGetDA(daapplication, i, &da); CHKERRQ(info);
    info = DAGetInfo(da,&dim,0,0,0,0,0,0,&dof,&s,0,&st); CHKERRQ(info);
    if (dim!=2){
      SETERRQ(1,"TAO DA ERROR: DA must have dimension of 2");}
    if (dof!=1){
      SETERRQ(1,"TAO DA ERROR: DA must have exactly 1 degree of freedom per node");}
    if (s!=1){
      SETERRQ(1,"TAO DA ERROR: DA stencil width must equal 1"); }
    if (st!=DA_STENCIL_BOX){
      SETERRQ(1,"TAO DA ERROR: DA stencil must be DA_STENCIL_BOX");}
  }
  PetscNew(TaoDA2D1DOFHCtx,&hctx);
  hctx->computeelementhessian=hess;
  hctx->elementhctx=ctx;
  hctx->elementhflops = flops;
  info = DAAppSetHessianRoutine(daapplication, TaoDA2dLoopHessian,  (void*)hctx);
  info = TaoAppSetDestroyRoutine(daapplication,TaoApplicationFreeMemory, (void*)hctx); CHKERRQ(info);
  info = PetscInfo(daapplication,"Set Hessian for TAO_APPLICATION object and create matrix.\n"); CHKERRQ(info);
  PetscFunctionReturn(0);
}

/*
typedef struct {
  / * Function Evaluation over entire single element  * /
  int  (*computeelementfunction)(int[2], double*,int, double*,void*); 
  int  (*computeelementgradient)(int[2], double*,int, double*,void*); 
  void *elementfctx;
  void *elementgctx;
  int elementfflops;
  int elementgflops;
} TaoPetscFDHessianAppCtx;
*/
/*
#undef __FUNCT__
#define __FUNCT__ "DAAppSetElementFunction"
int DAAppSetElementFunction(TAO_APPLICATION daapplication, int (*func)(int,int,double*,int,double*,void*), int flops, void *ctx){
  int info;
  PetscFunctionBegin;
  daapplication->computeelementfunction=func;
  daapplication->usrfctx=ctx;
  info = PetscLogInfo((daapplication->grid[0].da,"Set objective function pointer for TAO_APPLICATION object.\n")); CHKERRQ(info);
  PetscFunctionReturn(0);
} */



/*
#undef __FUNCT__
#define __FUNCT__ "DAAppSetElementGradient"

int DAAppSetElementGradient(TAO_APPLICATION daapplication, int (*grad)(int,int,double*,int,double*,void*), int flops, void *ctx){
  int info;
  PetscFunctionBegin;
  daapplication->computeelementgradient=grad;
  daapplication->usrgctx=ctx;
  info = PetscLogInfo((daapplication->grid[0].da,"Set gradient function pointer for TAO_APPLICATION object.\n")); CHKERRQ(info);
  PetscFunctionReturn(0);
}
*/


