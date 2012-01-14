#ifdef ad_GRAD_MAX
#undef ad_GRAD_MAX
#endif
#define ad_GRAD_MAX 4

#include "taodaapplication.h"   /* taodaapplication.h */
#include "taoapp.h"

//#include "ad_deriv.h"

/* Function Gradient Evaluation over single element  */
typedef struct {
  int (*computeadicfunctiongradient)(int[2],DERIV_TYPE[4],DERIV_TYPE*,void*);
  void *elementfgctx;
  int elementfgflops;
  DERIV_TYPE adX[4];
} TaoDA2D1DOFADICFGCtx;


#undef __FUNCT__
#define __FUNCT__ "TaoDA2dLoopADFunctionGradient"
static int TaoDA2dLoopADFunctionGradient(TAO_APPLICATION tao, DA da, Vec X, double *f, Vec G, void * ctx) {

  TaoDA2D1DOFADICFGCtx *myapp = (TaoDA2D1DOFADICFGCtx*) ctx;
  MPI_Comm comm;
  Vec localX, localG;
  int info, i, j, coor[2];
  int xs, xm, gxs, gxm, xe, ys, ym, gys, gym, ye;
  PetscScalar **x, **g;
  PetscScalar floc = 0.0;
  PetscScalar zero = 0.0;
  DERIV_TYPE adF,*adX=myapp->adX;

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

        DERIV_val(adX[0]) = x[j][i];
        DERIV_val(adX[1]) = x[j][i+1];
        DERIV_val(adX[2]) = x[j+1][i];
        DERIV_val(adX[3]) = x[j+1][i+1];
        coor[0] = i; coor[1] = j;

        info = myapp->computeadicfunctiongradient(coor,adX,&adF,myapp->elementfgctx);
	CHKERRQ(info);

        floc += DERIV_val(adF);

        g[j][i] += DERIV_grad(adF)[0];
        g[j][i+1] += DERIV_grad(adF)[1];
        g[j+1][i] += DERIV_grad(adF)[2];
        g[j+1][i+1] += DERIV_grad(adF)[3];
    }
  }

  PetscLogFlops((ye-ys)*(xe-xs)*(myapp->elementfgflops + 5));

  PetscObjectGetComm((PetscObject)da,&comm); CHKERRQ(info);
  info = MPI_Allreduce(&floc, f, 1, MPI_DOUBLE, MPI_SUM, comm); CHKERRQ(info);

  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecRestoreArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DALocalToGlobalBegin(da, localG, G); CHKERRQ(info);
  info = DALocalToGlobalEnd(da, localG, G); CHKERRQ(info);

  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localG); CHKERRQ(info);

  PetscFunctionReturn(0);
} /* TaoDA2dLoopADFunctionGradient */

#undef __FUNCT__
#define __FUNCT__ "TaoShutDownADIC"
static int TaoShutDownADICQQQQ(void* ctx){
  PetscFunctionBegin;
  ad_AD_Final();
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppSetADElementFunctionGradient"
/*@
   DAAppSetADElementFunctionGradient - Set routine that evaluates the
   local part of a function on a 2-dimensional DA with 1 degree of freedom. 

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapp - the TAO_DA_APPLICATION solver context
.  funcgrad - local function gradient routine
.  flops - the number of flops done performed in the funcgrad routine
-  fgctx - [optional] user-defined context for private data for the evaluation.

   Calling sequence of funcgrad:
$     int funcgrad(int coordinates[2], PetscScalar x[4], double *f, PetscScalar g[4], void* ctx)

+    coord - the global coordinates [i j] in each direction of the DA
.    x - the variables on the DA ( da[j][i], da[j][j+1], da[j+1][i], da[j+1][i+1] ) (bottom left, bottom right, top left, top right)
.    g - the ADIC differentiated objective function with respect to each variable
-    ctx - user defined context
   
   Note: This function requires ADIC to be installed and the ADIC-specific variables to be set in
         $TAO_DIR/bmake/packages.$PETSC_ARCH

   Level: intermediate

.keywords: DA, gradient, ADIC

.seealso:  DAAppSetObjectiveAndGradientRoutine(), DAAppSetElementObjectiveAndGradientRoutine()
@*/
int DAAppSetADElementFunctionGradient(TAO_APPLICATION daapplication, 
					 int (*funcgrad)(int[2],DERIV_TYPE[4],DERIV_TYPE*,void*), 
					 int flops, void *ctx){
  int i,n,info;
  int dim,dof,s;
  DAStencilType st;
  TaoDA2D1DOFADICFGCtx *fgctx;
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
  PetscNew(TaoDA2D1DOFADICFGCtx,&fgctx);
  //  ad_AD_Init(4);
  ad_AD_Init(ad_GRAD_MAX);
  ad_AD_SetIndepArray(fgctx->adX,4);
  ad_AD_SetIndepDone();
  fgctx->computeadicfunctiongradient=funcgrad;
  fgctx->elementfgctx=ctx;
  fgctx->elementfgflops=flops;
  info = DAAppSetObjectiveAndGradientRoutine(daapplication, TaoDA2dLoopADFunctionGradient, (void*)fgctx); 
  CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(daapplication,TaoApplicationFreeMemory, (void*)fgctx); CHKERRQ(info);
  info = TaoAppSetDestroyRoutine(daapplication, TaoShutDownADICQQQQ, 0); CHKERRQ(info);
  info = PetscInfo(daapplication,"Set objective function pointer for TAO_DA_APPLICATION object.\n"); CHKERRQ(info);
  PetscFunctionReturn(0);
}
