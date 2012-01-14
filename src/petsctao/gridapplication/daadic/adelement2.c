#ifdef ad_GRAD_MAX
#undef ad_GRAD_MAX
#endif
#define ad_GRAD_MAX 16

#include "taodaapplication.h"   /* taodaapplication.h */
#include "taoapp.h"

//#include "ad_deriv.h"

typedef struct {
  PetscScalar x,y,vpotx,vpoty;
} Field;


/* Function Gradient Evaluation over single element  */
typedef struct {
  int (*computeadicfunctiongradient)(int[2],DERIV_TYPE[16],DERIV_TYPE*,void*);
  void *elementfgctx;
  int elementfgflops;
  DERIV_TYPE adX[16];
} TaoDA2D4DOFADICFGCtx;


#undef __FUNCT__
#define __FUNCT__ "TaoDA2dLoopADFunctionGradient"
static int TaoDA2dLoopADFunctionGradient(TAO_APPLICATION tao, DA da, Vec X, double *f, Vec G, void * ctx) {

  TaoDA2D4DOFADICFGCtx *myapp = (TaoDA2D4DOFADICFGCtx*) ctx;
  MPI_Comm comm;
  Vec localX, localG;
  int info, i, j, coor[2];
  int xs, xm, mx, xe, ys, ym, ye, my;
  PetscScalar zero=0.0, floc = 0.0;
  DERIV_TYPE adF,*adX=myapp->adX;
  Field     **x,**g;

  info = DAGetLocalVector(da, &localX); CHKERRQ(info);
  info = DAGetLocalVector(da, &localG); CHKERRQ(info);
  info = VecSet(G, zero); CHKERRQ(info);
  info = VecSet(localG, zero); CHKERRQ(info);

  info = DAGetInfo(da,PETSC_NULL,&mx,&my,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,
		   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(info);

  info = DAGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(info);

  info = DAVecGetArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecGetArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DAGetCorners(da, &xs, &ys, PETSC_NULL, &xm, &ym, PETSC_NULL); CHKERRQ(info);

  xe = xs + xm;
  ye = ys + ym;

  for (j = ys; j < ye; j++) {
    for (i = xs; i < xe; i++) {

        DERIV_val(adX[0]) = x[j][i].x;
        DERIV_val(adX[1]) = x[j][i].y;
        DERIV_val(adX[2]) = x[j][i].vpotx;
        DERIV_val(adX[3]) = x[j][i].vpoty;

        DERIV_val(adX[4]) = x[j][i+1].x;
        DERIV_val(adX[5]) = x[j][i+1].y;
        DERIV_val(adX[6]) = x[j][i+1].vpotx;
        DERIV_val(adX[7]) = x[j][i+1].vpoty;

        DERIV_val(adX[8]) = x[j+1][i].x;
        DERIV_val(adX[9]) = x[j+1][i].y;
        DERIV_val(adX[10]) = x[j+1][i].vpotx;

        DERIV_val(adX[11]) = x[j+1][i].vpoty;

        DERIV_val(adX[12]) = x[j+1][i+1].x;
        DERIV_val(adX[13]) = x[j+1][i+1].y;
        DERIV_val(adX[14]) = x[j+1][i+1].vpotx;
        DERIV_val(adX[15]) = x[j+1][i+1].vpoty;

        coor[0] = i; coor[1] = j;

        info = myapp->computeadicfunctiongradient(coor,adX,&adF,myapp->elementfgctx);
	CHKERRQ(info);

        floc += DERIV_val(adF);

        g[j][i].x += DERIV_grad(adF)[0];
        g[j][i].y += DERIV_grad(adF)[1];
        g[j][i].vpotx += DERIV_grad(adF)[2];
        g[j][i].vpoty += DERIV_grad(adF)[3];

        g[j][i+1].x += DERIV_grad(adF)[4];
        g[j][i+1].y += DERIV_grad(adF)[5];
        g[j][i+1].vpotx += DERIV_grad(adF)[6];
        g[j][i+1].vpoty += DERIV_grad(adF)[7];

        g[j+1][i].x += DERIV_grad(adF)[8];
        g[j+1][i].y += DERIV_grad(adF)[9];
	g[j+1][i].vpotx += DERIV_grad(adF)[10];

	g[j+1][i].vpoty += DERIV_grad(adF)[11];

        g[j+1][i+1].x += DERIV_grad(adF)[12];
        g[j+1][i+1].y += DERIV_grad(adF)[13];
        g[j+1][i+1].vpotx += DERIV_grad(adF)[14];
        g[j+1][i+1].vpoty += DERIV_grad(adF)[15];

    }
  }


  info = DAVecRestoreArray(da, localX, (void**)&x); CHKERRQ(info);
  info = DAVecRestoreArray(da, localG, (void**)&g); CHKERRQ(info);

  info = DALocalToGlobalBegin(da, localG, G); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localX); CHKERRQ(info);
  info = DARestoreLocalVector(da, &localG); CHKERRQ(info);
  info = DALocalToGlobalEnd(da, localG, G); CHKERRQ(info);

  PetscLogFlops((ye-ys)*(xe-xs)*(myapp->elementfgflops + 33));
  PetscObjectGetComm((PetscObject)da,&comm); CHKERRQ(info);
  info = MPI_Allreduce(&floc, f, 1, MPI_DOUBLE, MPI_SUM, comm); CHKERRQ(info);

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
#define __FUNCT__ "DAAppSetADElementFunctionGradient2"
/*@C
   DAAppSetADElementFunctionGradient2 - Set routine that evaluates the
   local part of a function on a 2-dimensional periodic DA with 4 degrees of freedom. 

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapp - the TAO_DA_APPLICATION solver context
.  funcgrad - local function gradient routine
.  flops - the number of flops done performed in the funcgrad routine
-  fgctx - [optional] user-defined context for private data for the evaluation.

   Calling sequence of funcgrad:
$     int funcgrad(int coordinates[2], DERIV_TYPE x[16], DERIV_TYPE *f, void* ctx)

+    coord - the global coordinates [i j] in each direction of the DA
.    x - the variables on the DA ( da[j][i], da[j][j+1], da[j+1][i], da[j+1][i+1] ) (bottom left, bottom right, top left, top right)
.    g - the ADIC differentiated objective function with respect to each variable
-    ctx - user defined context

   Calling sequence of func before calling ADIC
$     int func(int coordinates[2], double x[16], double *f, void* ctx)

+    coord - the global coordinates [i j] in each direction of the DA
.    x - the variables on the DA ( da[j][i], da[j][j+1], da[j+1][i], da[j+1][i+1] ) (bottom left, bottom right, top left, top right)
.    f - the objective function
-    ctx - user defined context


   Note: This function requires ADIC to be installed and the ADIC-specific variables to be set in
         $TAO_DIR/bmake/packages.$PETSC_ARCH

   Level: intermediate

Concepts: DA, gradient, ADIC

.seealso: DAAppSetObjectiveAndGradientRoutine(), DAAppSetADElementFunctionGradient()
@*/
int DAAppSetADElementFunctionGradient2(TAO_APPLICATION daapplication, 
					 int (*funcgrad)(int[2],DERIV_TYPE[16],DERIV_TYPE*,void*), 
					 int flops, void *ctx){
  int i,n,info;
  int dim,dof,s;
  DAStencilType st;
  TaoDA2D4DOFADICFGCtx *fgctx;
  DA da;

  PetscFunctionBegin;
  info=DAAppGetNumberOfDAGrids(daapplication,&n); CHKERRQ(info);
  for (i=0;i<n;i++){
    info = DAAppGetDA(daapplication, i, &da); CHKERRQ(info);
    info = DAGetInfo(da,&dim,0,0,0,0,0,0,&dof,&s,0,&st); CHKERRQ(info);
    if (dim!=2){
      SETERRQ(1,"TAO DA ERROR: DA must have dimension of 2");}
    if (dof!=4){
      SETERRQ(1,"TAO DA ERROR: DA must have exactly 4 degrees of freedom per node");}
    if (s!=1){
      SETERRQ(1,"TAO DA ERROR: DA stencil width must equal 1"); }
    //    if (st!=DA_STENCIL_BOX){
    //      SETERRQ(1,"TAO DA ERROR: DA stencil must be DA_STENCIL_BOX");}
  }
  PetscNew(TaoDA2D4DOFADICFGCtx,&fgctx);
  ad_AD_Init(16);
  //  ad_AD_Init(ad_GRAD_MAX);
  ad_AD_SetIndepArray(fgctx->adX,16);
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

