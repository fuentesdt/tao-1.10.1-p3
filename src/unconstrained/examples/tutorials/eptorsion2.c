/*$Id$*/

/* Program usage: mpirun -np <proc> eptorsion2 [-help] [all TAO options] */

/* ----------------------------------------------------------------------

  Elastic-plastic torsion problem.  

  The elastic plastic torsion problem arises from the determination 
  of the stress field on an infinitely long cylindrical bar, which is
  equivalent to the solution of the following problem:

  min{ .5 * integral(||gradient(v(x))||^2 dx) - C * integral(v(x) dx)}
  
  where C is the torsion angle per unit length.

  The uniprocessor version of this code is eptorsion1.c; the Fortran 
  version of this code is eptorsion2f.F.

  This application solves the problem without calculating hessians 
---------------------------------------------------------------------- */

/*
  Include "tao.h" so that we can use TAO solvers.  Note that this 
  file automatically includes files for lower-level support, such as those
  provided by the PETSc library:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - sysem routines        petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
  Include "petscda.h" so that we can use distributed arrays (DAs) for managing
  the parallel mesh.
*/

#include "tao.h"
#include "petscda.h"

static  char help[] = 
"Demonstrates use of the TAO package to solve \n\
unconstrained minimization problems in parallel.  This example is based on \n\
the Elastic-Plastic Torsion (dept) problem from the MINPACK-2 test suite.\n\
The command line options are:\n\
  -mx <xg>, where <xg> = number of grid points in the 1st coordinate direction\n\
  -my <yg>, where <yg> = number of grid points in the 2nd coordinate direction\n\
  -par <param>, where <param> = angle of twist per unit length\n\n";

/*T
   Concepts: TAO - Solving an unconstrained minimization problem
   Routines: TaoInitialize(); TaoFinalize();              
   Routines: TaoApplicationCreate(); TaoAppDestroy();
   Routines: TaoCreate(); TaoDestroy(); 
   Routines: TaoAppSetObjectiveAndGradientRoutine();
   Routines: TaoAppSetHessianMat(); TaoAppSetHessianRoutine();
   Routines: TaoSetOptions();
   Routines: TaoAppSetInitialSolutionVec();
   Routines: TaoSolve(); TaoView(); TaoAppGetKSP();
   Routines: TaoGetSolutionStatus(); 
   Processors: n
T*/



/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines, FormFunction() and
   FormGradient().
*/
typedef struct {
  /* parameters */
   PetscInt           mx, my;         /* global discretization in x- and y-directions */
   PetscReal        param;          /* nonlinearity parameter */

  /* work space */
   Vec           localX;         /* local vectors */
   DA            da;             /* distributed array data structure */
} AppCtx;

/* -------- User-defined Routines --------- */

int FormInitialGuess(AppCtx*,Vec);
int FormFunctionGradient(TAO_APPLICATION,Vec,double*,Vec,void*);
int ComputeHessian(TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure*,void*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  int             info;                  /* used to check for functions returning nonzeros */
  Vec             x;                     /* solution, gradient vectors */
  Mat             H;                     /* Hessian matrix */
  PetscInt        Nx, Ny;                /* number of processes in x- and y-directions */
  TAO_SOLVER      tao;                   /* TAO_SOLVER solver context */
  TAO_APPLICATION torsionapp;            /* TAO application context */
  TaoTerminateReason reason;
  PetscTruth      flg;
  PetscInt        iter;                  /* iteration information */
  double          ff,gnorm;
  AppCtx          user;                  /* application context */
  KSP             ksp;

  /* Initialize TAO, PETSc contexts */
  info = PetscInitialize(&argc,&argv,(char *)0,help);
  info = TaoInitialize(&argc,&argv,(char *)0,help);

  /* Specify default parameters */
  user.param = 5.0; user.mx = 10; user.my = 10;
  Nx = Ny = PETSC_DECIDE;

  /* Check for any command line arguments that override defaults */
  info = PetscOptionsGetReal(TAO_NULL,"-par",&user.param,&flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL,"-my",&user.my,&flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL,"-mx",&user.mx,&flg); CHKERRQ(info);

  PetscPrintf(PETSC_COMM_WORLD,"\n---- Elastic-Plastic Torsion Problem -----\n");
  PetscPrintf(PETSC_COMM_WORLD,"mx: %d     my: %d   \n\n",user.mx,user.my);  

  /* Set up distributed array */
  info = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
		    user.mx,user.my,Nx,Ny,1,1,TAO_NULL,TAO_NULL,
		    &user.da); CHKERRQ(info);

  /* Create vectors */
  info = DACreateGlobalVector(user.da,&x); CHKERRQ(info);

  info = DACreateLocalVector(user.da,&user.localX); CHKERRQ(info);

  /* Create Hessian */
  info = DAGetMatrix(user.da,MATAIJ,&H); CHKERRQ(info);
  info = MatSetOption(H,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(info);

  /* The TAO code begins here */

  /* Create TAO solver and set desired solution method */
  info = TaoCreate(MPI_COMM_WORLD,"tao_cg",&tao); CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_WORLD,&torsionapp); CHKERRQ(info);

  /* Set initial solution guess */
  info = FormInitialGuess(&user,x); CHKERRQ(info);
  info = TaoAppSetInitialSolutionVec(torsionapp,x); CHKERRQ(info);

  /* Set routine for function and gradient evaluation */
  info = TaoAppSetObjectiveAndGradientRoutine(torsionapp,FormFunctionGradient,(void *)&user); CHKERRQ(info);

  info = TaoAppSetHessianMat(torsionapp,H,H); CHKERRQ(info);
  info = TaoAppSetHessianRoutine(torsionapp,ComputeHessian,(void*)&user); CHKERRQ(info);

  info = TaoAppGetKSP(torsionapp,&ksp); CHKERRQ(info);
  if (ksp) {                                              /* Modify the PETSc KSP structure */
    info = KSPSetType(ksp,KSPCG); CHKERRQ(info);
  }

  /* Check for any TAO command line options */
  info = TaoSetOptions(torsionapp,tao); CHKERRQ(info);

  /* SOLVE THE APPLICATION */
  info = TaoSolveApplication(torsionapp,tao);  CHKERRQ(info);

  /* Get information on termination */
  info = TaoGetSolutionStatus(tao,&iter,&ff,&gnorm,0,0,&reason); CHKERRQ(info);
  if (reason <= 0){
    PetscPrintf(MPI_COMM_WORLD, "Try another method! Iterations: %d, f: %4.2e, residual: %4.2e\n",
		iter,ff,gnorm); CHKERRQ(info); 
  }  

  /* 
     To View TAO solver information use
      info = TaoView(tao); CHKERRQ(info);
  */

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(torsionapp);  CHKERRQ(info);

  /* Free PETSc data structures */
  info = VecDestroy(x); CHKERRQ(info);
  info = MatDestroy(H); CHKERRQ(info);

  info = VecDestroy(user.localX); CHKERRQ(info);
  info = DADestroy(user.da); CHKERRQ(info);


  /* Finalize TAO and PETSc*/
  TaoFinalize();
  PetscFinalize();

  return 0;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
/*
    FormInitialGuess - Computes an initial approximation to the solution.

    Input Parameters:
.   user - user-defined application context
.   X    - vector

    Output Parameters:
    X    - vector
*/
int FormInitialGuess(AppCtx *user,Vec X)
{
  int    info;
  PetscInt   i, j, k, mx = user->mx, my = user->my;
  PetscInt   xs, ys, xm, ym, gxm, gym, gxs, gys, xe, ye;
  PetscReal hx = 1.0/(mx+1), hy = 1.0/(my+1), temp, val;

  /* Get local mesh boundaries */
  info = DAGetCorners(user->da,&xs,&ys,TAO_NULL,&xm,&ym,TAO_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(user->da,&gxs,&gys,TAO_NULL,&gxm,&gym,TAO_NULL); CHKERRQ(info);

  /* Compute initial guess over locally owned part of mesh */
  xe = xs+xm;
  ye = ys+ym;
  for (j=ys; j<ye; j++) {  /*  for (j=0; j<my; j++) */
    temp = TaoMin(j+1,my-j)*hy;
    for (i=xs; i<xe; i++) {  /*  for (i=0; i<mx; i++) */
      k   = (j-gys)*gxm + i-gxs;
      val = TaoMin((TaoMin(i+1,mx-i))*hx,temp);
      info = VecSetValuesLocal(X,1,&k,&val,ADD_VALUES); CHKERRQ(info);
    }
  }
  info = VecAssemblyBegin(X); CHKERRQ(info);
  info = VecAssemblyEnd(X); CHKERRQ(info);

  return 0;
}


/* ------------------------------------------------------------------ */
#undef __FUNCT__
#define __FUNCT__ "FormFunctionGradient"
/* 
   FormFunctionGradient - Evaluates the function and corresponding gradient.
    
   Input Parameters:
   taoapp - the TAO_APPLICATION context
   X   - the input vector 
   ptr - optional user-defined context, as set by TaoSetFunction()

   Output Parameters:
   f   - the newly evaluated function
   G   - the newly evaluated gradient
*/
int FormFunctionGradient(TAO_APPLICATION taoapp,Vec X,double *f,Vec G,void *ptr){

  AppCtx *user = (AppCtx *)ptr;
  int    info;
  PetscInt i,j,k,ind;
  PetscInt xe,ye,xsm,ysm,xep,yep;
  PetscInt xs, ys, xm, ym, gxm, gym, gxs, gys;
  PetscInt mx = user->mx, my = user->my;
  PetscReal three = 3.0, zero = 0.0, *x, floc, cdiv3 = user->param/three;
  PetscReal p5 = 0.5, area, val, flin, fquad;
  PetscReal v,vb,vl,vr,vt,dvdx,dvdy;
  PetscReal hx = 1.0/(user->mx + 1); 
  PetscReal hy = 1.0/(user->my + 1);  
  Vec    localX = user->localX;


  /* Initialize */
  flin = fquad = zero;

  info = VecSet(G, zero); CHKERRQ(info);
  /*
     Scatter ghost points to local vector,using the 2-step process
        DAGlobalToLocalBegin(),DAGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  info = DAGlobalToLocalBegin(user->da,X,INSERT_VALUES,localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(user->da,X,INSERT_VALUES,localX); CHKERRQ(info);

  /* Get pointer to vector data */
  info = VecGetArray(localX,&x); CHKERRQ(info);

  /* Get local mesh boundaries */
  info = DAGetCorners(user->da,&xs,&ys,TAO_NULL,&xm,&ym,TAO_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(user->da,&gxs,&gys,TAO_NULL,&gxm,&gym,TAO_NULL); CHKERRQ(info);

  /* Set local loop dimensions */
  xe = xs+xm;
  ye = ys+ym;
  if (xs == 0)  xsm = xs-1;
  else          xsm = xs;
  if (ys == 0)  ysm = ys-1;
  else          ysm = ys;
  if (xe == mx) xep = xe+1;
  else          xep = xe;
  if (ye == my) yep = ye+1;
  else          yep = ye;

  /* Compute local gradient contributions over the lower triangular elements */
  for (j=ysm; j<ye; j++) {  /*  for (j=-1; j<my; j++) */
    for (i=xsm; i<xe; i++) {  /*   for (i=-1; i<mx; i++) */
      k = (j-gys)*gxm + i-gxs;
      v = zero;
      vr = zero;
      vt = zero;
      if (i >= 0 && j >= 0) v = x[k];
      if (i < mx-1 && j > -1) vr = x[k+1];
      if (i > -1 && j < my-1) vt = x[k+gxm];
      dvdx = (vr-v)/hx;
      dvdy = (vt-v)/hy;
      if (i != -1 && j != -1) {
        ind = k; val = - dvdx/hx - dvdy/hy - cdiv3;
        info = VecSetValuesLocal(G,1,&k,&val,ADD_VALUES); CHKERRQ(info);
      }
      if (i != mx-1 && j != -1) {
        ind = k+1; val =  dvdx/hx - cdiv3;
        info = VecSetValuesLocal(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      if (i != -1 && j != my-1) {
        ind = k+gxm; val = dvdy/hy - cdiv3;
        info = VecSetValuesLocal(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      fquad += dvdx*dvdx + dvdy*dvdy;
      flin -= cdiv3 * (v + vr + vt);
    }
  }

  /* Compute local gradient contributions over the upper triangular elements */
  for (j=ys; j<yep; j++) { /*  for (j=0; j<=my; j++) */
    for (i=xs; i<xep; i++) {  /*   for (i=0; i<=mx; i++) */
      k = (j-gys)*gxm + i-gxs;
      vb = zero;
      vl = zero;
      v  = zero;
      if (i < mx && j > 0) vb = x[k-gxm];
      if (i > 0 && j < my) vl = x[k-1];
      if (i < mx && j < my) v = x[k];
      dvdx = (v-vl)/hx;
      dvdy = (v-vb)/hy;
      if (i != mx && j != 0) {
        ind = k-gxm; val = - dvdy/hy - cdiv3;
        info = VecSetValuesLocal(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      if (i != 0 && j != my) {
        ind = k-1; val =  - dvdx/hx - cdiv3;
        info = VecSetValuesLocal(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      if (i != mx && j != my) {
        ind = k; val =  dvdx/hx + dvdy/hy - cdiv3;
        info = VecSetValuesLocal(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      fquad += dvdx*dvdx + dvdy*dvdy;
      flin -= cdiv3 * (vb + vl + v);
    }
  }


  /* Restore vector */
  info = VecRestoreArray(localX,&x); CHKERRQ(info);

  /* Assemble gradient vector */
  info = VecAssemblyBegin(G); CHKERRQ(info);
  info = VecAssemblyEnd(G); CHKERRQ(info);

  /* Scale the gradient */
  area = p5*hx*hy;
  floc = area * (p5 * fquad + flin);
  info = VecScale(G, area); CHKERRQ(info);

  /* Sum function contributions from all processes */
  MPI_Allreduce((void*)&floc,(void*)f,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  info=PetscLogFlops((ye-ysm)*(xe-xsm)*20+(xep-xs)*(yep-ys)*16); CHKERRQ(info);
  
  return 0;
}



#undef __FUNCT__
#define __FUNCT__ "ComputeHessian"
int ComputeHessian(TAO_APPLICATION taoapp, Vec X, Mat *H, Mat *Hpre, MatStructure *flag, void*ctx){

  AppCtx *user= (AppCtx*) ctx;
  int info;
  PetscInt i,j,k;
  PetscInt col[5],row;
  PetscInt xs,xm,gxs,gxm,ys,ym,gys,gym;
  PetscReal v[5];
  PetscReal hx=1.0/(user->mx+1), hy=1.0/(user->my+1), hxhx=1.0/(hx*hx), hyhy=1.0/(hy*hy), area=0.5*hx*hy;
  Mat A=*H;

  /* Compute the quadratic term in the objective function */  

  /*
     Get local grid boundaries
  */

  info = DAGetCorners(user->da,&xs,&ys,TAO_NULL,&xm,&ym,TAO_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(user->da,&gxs,&gys,TAO_NULL,&gxm,&gym,TAO_NULL); CHKERRQ(info);

  for (j=ys; j<ys+ym; j++){
    
    for (i=xs; i< xs+xm; i++){

      row=(j-gys)*gxm + (i-gxs);

      k=0;
      if (j>gys){ 
        v[k]=-2*hyhy; col[k]=row - gxm; k++;
      }

      if (i>gxs){
        v[k]= -2*hxhx; col[k]=row - 1; k++;
      }

      v[k]= 4.0*(hxhx+hyhy); col[k]=row; k++;

      if (i+1 < gxs+gxm){
        v[k]= -2.0*hxhx; col[k]=row+1; k++;
      }

      if (j+1 <gys+gym){
        v[k]= -2*hyhy; col[k] = row+gxm; k++;
      }

      info = MatSetValuesLocal(A,1,&row,k,col,v,INSERT_VALUES); CHKERRQ(info);

    }
  }
  /* 
     Assemble matrix, using the 2-step process:
     MatAssemblyBegin(), MatAssemblyEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  info = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  /*
    Tell the matrix we will never add a new nonzero location to the
    matrix. If we do it will generate an error.
  */
  info = MatScale(A,area); CHKERRQ(info);
  info = MatSetOption(A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(info);
  info = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(info);

  info = PetscLogFlops(9*xm*ym+49*xm); CHKERRQ(info);

  return 0;
}













