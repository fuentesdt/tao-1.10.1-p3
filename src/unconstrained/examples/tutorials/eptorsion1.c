/*$Id$*/

/* Program usage: mpirun -np 1 eptorsion1 [-help] [all TAO options] */

/* ----------------------------------------------------------------------

  Elastic-plastic torsion problem.  

  The elastic plastic torsion problem arises from the determination 
  of the stress field on an infinitely long cylindrical bar, which is
  equivalent to the solution of the following problem:

  min{ .5 * integral(||gradient(v(x))||^2 dx) - C * integral(v(x) dx)}
  
  where C is the torsion angle per unit length.

  The multiprocessor version of this code is eptorsion2.c.

---------------------------------------------------------------------- */

/*
  Include "tao.h" so that we can use TAO solvers.  Note that this 
  file automatically includes files for lower-level support, such as those
  provided by the PETSc library:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - sysem routines        petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/

#include "tao.h"


static  char help[]=
"Demonstrates use of the TAO package to solve \n\
unconstrained minimization problems on a single processor.  This example \n\
is based on the Elastic-Plastic Torsion (dept) problem from the MINPACK-2 \n\
test suite.\n\
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
   Routines: TaoSolveApplication();
   Routines: TaoGetSolutionStatus(); TaoAppGetKSP();
   Processors: 1
T*/ 

/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines, FormFunction(),
   FormGradient(), and FormHessian().
*/

typedef struct {
   PetscReal  param;      /* nonlinearity parameter */
   PetscInt   mx, my;     /* discretization in x- and y-directions */
   PetscInt   ndim;       /* problem dimension */
   Vec     s, y, xvec; /* work space for computing Hessian */
   PetscReal  hx, hy;     /* mesh spacing in x- and y-directions */
} AppCtx;

/* -------- User-defined Routines --------- */

int FormInitialGuess(AppCtx*,Vec);
int FormFunction(TAO_APPLICATION,Vec,double*,void*);
int FormGradient(TAO_APPLICATION,Vec,Vec,void*);
int FormHessian(TAO_APPLICATION,Vec,Mat*,Mat*, MatStructure *,void*);
int HessianProductMat(Mat,Vec,Vec);
int HessianProduct(void*,Vec,Vec);
int MatrixFreeHessian(TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure*,void*);
int FormFunctionGradient(TAO_APPLICATION,Vec,double *,Vec,void *);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  int         info;                /* used to check for functions returning nonzeros */
  PetscInt      mx=10;               /* discretization in x-direction */
  PetscInt      my=10;               /* discretization in y-direction */
  Vec         x;                   /* solution, gradient vectors */
  PetscTruth  flg;                 /* A return value when checking for use options */
  TAO_SOLVER  tao;                 /* TAO_SOLVER solver context */
  TAO_APPLICATION eptorsionapp;    /* The PETSc application */
  Mat         H;                   /* Hessian matrix */
  TaoInt      iter;                /* iteration information */
  double      ff,gnorm;
  TaoTerminateReason reason;        
  KSP         ksp;                 /* PETSc Krylov subspace solver */
  PC          pc;                  /* PETSc preconditioner */
  AppCtx      user;                /* application context */
  int         size;                /* number of processes */
  double      one=1.0;


  /* Initialize TAO,PETSc */
  PetscInitialize(&argc,&argv,(char *)0,help);
  TaoInitialize(&argc,&argv,(char *)0,help);

  /* Optional:  Check  that only one processor is being used. */
  info = MPI_Comm_size(MPI_COMM_WORLD,&size); CHKERRQ(info);
  if (size >1) {
    PetscPrintf(PETSC_COMM_SELF,"This example is intended for single processor use!\n");
    PetscPrintf(PETSC_COMM_SELF,"Try the example eptorsion2!\n");
    SETERRQ(1,"Incorrect number of processors");
  }

  /* Specify default parameters for the problem, check for command-line overrides */
  user.param = 5.0;
  info = PetscOptionsGetInt(TAO_NULL,"-my",&my,&flg); CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL,"-mx",&mx,&flg); CHKERRQ(info);
  info = PetscOptionsGetReal(TAO_NULL,"-par",&user.param,&flg); CHKERRQ(info);


  PetscPrintf(PETSC_COMM_SELF,"\n---- Elastic-Plastic Torsion Problem -----\n");
  PetscPrintf(PETSC_COMM_SELF,"mx: %d     my: %d   \n\n",mx,my);  
  user.ndim = mx * my; user.mx = mx; user.my = my;

  user.hx = one/(mx+1); user.hy = one/(my+1);


  /* Allocate vectors */
  info = VecCreateSeq(PETSC_COMM_SELF,user.ndim,&user.y); CHKERRQ(info);
  info = VecDuplicate(user.y,&user.s); CHKERRQ(info);
  info = VecDuplicate(user.y,&x); CHKERRQ(info);

  /* The TAO code begins here */

  /* Create TAO solver and set desired solution method */
  info = TaoCreate(PETSC_COMM_SELF,"tao_lmvm",&tao); CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_SELF,&eptorsionapp); CHKERRQ(info);

  /* Set solution vector with an initial guess */
  info = FormInitialGuess(&user,x); CHKERRQ(info);
  info = TaoAppSetInitialSolutionVec(eptorsionapp,x); CHKERRQ(info);

  /* Set routine for function and gradient evaluation */
  info = TaoAppSetObjectiveAndGradientRoutine(eptorsionapp,FormFunctionGradient,(void *)&user); CHKERRQ(info);

  /* From command line options, determine if using matrix-free hessian */
  info = PetscOptionsHasName(TAO_NULL,"-my_tao_mf",&flg); CHKERRQ(info);
  if (flg) {
    info = MatCreateShell(PETSC_COMM_SELF,user.ndim,user.ndim,user.ndim,
                          user.ndim,(void*)&user,&H); CHKERRQ(info);
    info = MatShellSetOperation(H,MATOP_MULT,(void(*)())HessianProductMat); CHKERRQ
(info);
    info = MatSetOption(H,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(info);

    info = TaoAppSetHessianRoutine(eptorsionapp,MatrixFreeHessian,(void *)&user); CHKERRQ(info);
    info = TaoAppSetHessianMat(eptorsionapp,H,H); CHKERRQ(info);

    /* Set null preconditioner.  Alternatively, set user-provided 
       preconditioner or explicitly form preconditioning matrix */
    info = TaoAppGetKSP(eptorsionapp,&ksp); CHKERRQ(info);
    if (ksp){
      info = KSPGetPC(ksp,&pc); CHKERRQ(info);
      info = PCSetType(pc,PCNONE); CHKERRQ(info);
    }
  } else {

    info = MatCreateSeqAIJ(PETSC_COMM_SELF,user.ndim,user.ndim,5,TAO_NULL,&H); CHKERRQ(info);
    info = MatSetOption(H,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(info);

    info = TaoAppSetHessianRoutine(eptorsionapp,FormHessian,(void *)&user); CHKERRQ(info);
    info = TaoAppSetHessianMat(eptorsionapp,H,H); CHKERRQ(info); 
  }

  /* Check for any TAO command line options */
  info = TaoSetOptions(eptorsionapp,tao); CHKERRQ(info);


  /* Modify the PETSc KSP structure */
  info = TaoAppGetKSP(eptorsionapp,&ksp); CHKERRQ(info);
  if (ksp) {                                              
    info = KSPSetType(ksp,KSPCG); CHKERRQ(info);
  }

  /* SOLVE THE APPLICATION */
  info = TaoSolveApplication(eptorsionapp,tao);  CHKERRQ(info);
  if (ksp) {
    KSPView(ksp,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(info);
  }

  /* 
     To View TAO solver information use
      info = TaoView(tao); CHKERRQ(info);
  */

  /* Get information on termination */
  info = TaoGetSolutionStatus(tao,&iter,&ff,&gnorm,0,0,&reason); CHKERRQ(info);
  if (reason <= 0){
    PetscPrintf(PETSC_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
    PetscPrintf(PETSC_COMM_WORLD,"Iter: %d,   f: %4.2e,  residual: %4.2e\n",iter,ff,gnorm); CHKERRQ(info);
  }

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(eptorsionapp); CHKERRQ(info);

  /* Free PETSc data structures */
  info = VecDestroy(user.s); CHKERRQ(info);
  info = VecDestroy(user.y); CHKERRQ(info);
  info = VecDestroy(x); CHKERRQ(info);
  info = MatDestroy(H); CHKERRQ(info);

  /* Finalize TAO, PETSc */
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
.   X    - vector
*/
int FormInitialGuess(AppCtx *user,Vec X)
{
  double hx = user->hx, hy = user->hy, temp;
  PetscScalar val;
  int    info;
  PetscInt i, j, k, nx = user->mx, ny = user->my;

  /* Compute initial guess */
  for (j=0; j<ny; j++) {
    temp = TaoMin(j+1,ny-j)*hy;
    for (i=0; i<nx; i++) {
      k   = nx*j + i;
      val = TaoMin((TaoMin(i+1,nx-i))*hx,temp);
      info = VecSetValues(X,1,&k,&val,ADD_VALUES); CHKERRQ(info);
    }
  }
  info = VecAssemblyBegin(X); CHKERRQ(info);
  info = VecAssemblyEnd(X); CHKERRQ(info);
  return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunctionGradient"
/* 
   FormFunctionGradient - Evaluates the function and corresponding gradient.
    
   Input Parameters:
   tao - the TAO_APPLICATION context
   X   - the input vector 
   ptr - optional user-defined context, as set by TaoSetFunction()

   Output Parameters:
   f   - the newly evaluated function
   G   - the newly evaluated gradient
*/
int FormFunctionGradient(TAO_APPLICATION tao,Vec X,double *f,Vec G,void *ptr)
{
  int info;
  info = FormFunction(tao,X,f,ptr);CHKERRQ(info);
  info = FormGradient(tao,X,G,ptr);CHKERRQ(info);
  return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
/* 
   FormFunction - Evaluates the function, f(X).

   Input Parameters:
.  taoapp - the TAO_APPLICATION context
.  X   - the input vector 
.  ptr - optional user-defined context, as set by TaoSetFunction()

   Output Parameters:
.  f    - the newly evaluated function
*/
int FormFunction(TAO_APPLICATION taoapp,Vec X,double *f,void *ptr)
{
  AppCtx *user = (AppCtx *) ptr;
  double hx = user->hx, hy = user->hy, area, three = 3.0, p5 = 0.5;
  double zero = 0.0, vb, vl, vr, vt, dvdx, dvdy, flin = 0.0, fquad = 0.0;
  double v, cdiv3 = user->param/three;
  PetscScalar *x;
  int    info;
  PetscInt  nx = user->mx, ny = user->my, i, j, k;

  /* Get pointer to vector data */
  info = VecGetArray(X,&x); CHKERRQ(info);

  /* Compute function contributions over the lower triangular elements */
  for (j=-1; j<ny; j++) {
    for (i=-1; i<nx; i++) {
      k = nx*j + i;
      v = zero;
      vr = zero;
      vt = zero;
      if (i >= 0 && j >= 0) v = x[k];
      if (i < nx-1 && j > -1) vr = x[k+1];
      if (i > -1 && j < ny-1) vt = x[k+nx];
      dvdx = (vr-v)/hx;
      dvdy = (vt-v)/hy;
      fquad += dvdx*dvdx + dvdy*dvdy;
      flin -= cdiv3*(v+vr+vt);
    }
  }

  /* Compute function contributions over the upper triangular elements */
  for (j=0; j<=ny; j++) {
    for (i=0; i<=nx; i++) {
      k = nx*j + i;
      vb = zero;
      vl = zero;
      v = zero;
      if (i < nx && j > 0) vb = x[k-nx];
      if (i > 0 && j < ny) vl = x[k-1];
      if (i < nx && j < ny) v = x[k];
      dvdx = (v-vl)/hx;
      dvdy = (v-vb)/hy;
      fquad = fquad + dvdx*dvdx + dvdy*dvdy;
      flin = flin - cdiv3*(vb+vl+v);
    }
  }

  /* Restore vector */
  info = VecRestoreArray(X,&x); CHKERRQ(info);

  /* Scale the function */
  area = p5*hx*hy;
  *f = area*(p5*fquad+flin);

  info = PetscLogFlops(nx*ny*24); CHKERRQ(info);
  return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormGradient"
/*  
    FormGradient - Evaluates the gradient, G(X).              

    Input Parameters:
.   taoapp  - the TAO_APPLICATION context
.   X    - input vector
.   ptr  - optional user-defined context
    
    Output Parameters:
.   G - vector containing the newly evaluated gradient
*/
int FormGradient(TAO_APPLICATION taoapp,Vec X,Vec G,void *ptr)
{
  AppCtx *user = (AppCtx *) ptr;
  PetscScalar zero=0.0, p5=0.5, three = 3.0, area, val, *x;
  int    info;
  PetscInt nx = user->mx, ny = user->my, ind, i, j, k;
  double hx = user->hx, hy = user->hy;
  double vb, vl, vr, vt, dvdx, dvdy;
  double v, cdiv3 = user->param/three;

  /* Initialize gradient to zero */
  info = VecSet(G, zero); CHKERRQ(info);

  /* Get pointer to vector data */
  info = VecGetArray(X,&x); CHKERRQ(info);

  /* Compute gradient contributions over the lower triangular elements */
  for (j=-1; j<ny; j++) {
    for (i=-1; i<nx; i++) {
      k  = nx*j + i;
      v  = zero;
      vr = zero;
      vt = zero;
      if (i >= 0 && j >= 0)    v = x[k];
      if (i < nx-1 && j > -1) vr = x[k+1];
      if (i > -1 && j < ny-1) vt = x[k+nx];
      dvdx = (vr-v)/hx;
      dvdy = (vt-v)/hy;
      if (i != -1 && j != -1) {
        ind = k; val = - dvdx/hx - dvdy/hy - cdiv3;
        info = VecSetValues(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      if (i != nx-1 && j != -1) {
        ind = k+1; val =  dvdx/hx - cdiv3;
        info = VecSetValues(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      if (i != -1 && j != ny-1) {
        ind = k+nx; val = dvdy/hy - cdiv3;
        info = VecSetValues(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
    }
  }

  /* Compute gradient contributions over the upper triangular elements */
  for (j=0; j<=ny; j++) {
    for (i=0; i<=nx; i++) {
      k = nx*j + i;
      vb = zero;
      vl = zero;
      v  = zero;
      if (i < nx && j > 0) vb = x[k-nx];
      if (i > 0 && j < ny) vl = x[k-1];
      if (i < nx && j < ny) v = x[k];
      dvdx = (v-vl)/hx;
      dvdy = (v-vb)/hy;
      if (i != nx && j != 0) {
        ind = k-nx; val = - dvdy/hy - cdiv3;
        info = VecSetValues(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      if (i != 0 && j != ny) {
        ind = k-1; val =  - dvdx/hx - cdiv3;
        info = VecSetValues(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
      if (i != nx && j != ny) {
        ind = k; val =  dvdx/hx + dvdy/hy - cdiv3;
        info = VecSetValues(G,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
      }
    }
  }

  /* Restore vector */
  info = VecRestoreArray(X,&x); CHKERRQ(info);

  /* Assemble gradient vector */
  info = VecAssemblyBegin(G); CHKERRQ(info);
  info = VecAssemblyEnd(G); CHKERRQ(info);

  /* Scale the gradient */
  area = p5*hx*hy;
  info = VecScale(G, area); CHKERRQ(info);
  
  info = PetscLogFlops(nx*ny*24); CHKERRQ(info);
  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormHessian"
/* 
   FormHessian - Forms the Hessian matrix.

   Input Parameters:
.  taoapp - the TAO_APPLICATION context
.  X    - the input vector
.  ptr  - optional user-defined context, as set by TaoSetHessian()
   
   Output Parameters:
.  H     - Hessian matrix
.  PrecH - optionally different preconditioning Hessian
.  flag  - flag indicating matrix structure

   Notes:
   This routine is intended simply as an example of the interface
   to a Hessian evaluation routine.  Since this example compute the
   Hessian a column at a time, it is not particularly efficient and
   is not recommended.
*/
int FormHessian(TAO_APPLICATION taoapp,Vec X,Mat *HH,Mat *Hpre, MatStructure *flg, void *ptr)
{
  AppCtx     *user = (AppCtx *) ptr;
  int info;
  PetscInt   i,j, ndim = user->ndim;
  PetscScalar  *y, zero = 0.0, one = 1.0;
  Mat H=*HH;
  *Hpre = H;
  PetscTruth assembled;

  /* Set location of vector */
  user->xvec = X;

  /* Initialize Hessian entries and work vector to zero */
  info = MatAssembled(H,&assembled); CHKERRQ(info);
  if (assembled){info = MatZeroEntries(H);  CHKERRQ(info);}

  info = VecSet(user->s, zero); CHKERRQ(info);

  /* Loop over matrix columns to compute entries of the Hessian */
  for (j=0; j<ndim; j++) {

    info = VecSetValues(user->s,1,&j,&one,INSERT_VALUES); CHKERRQ(info);
    info = VecAssemblyBegin(user->s); CHKERRQ(info);
    info = VecAssemblyEnd(user->s); CHKERRQ(info);

    info = HessianProduct(ptr,user->s,user->y); CHKERRQ(info);

    info = VecSetValues(user->s,1,&j,&zero,INSERT_VALUES); CHKERRQ(info);
    info = VecAssemblyBegin(user->s); CHKERRQ(info);
    info = VecAssemblyEnd(user->s); CHKERRQ(info);

    info = VecGetArray(user->y,&y); CHKERRQ(info);
    for (i=0; i<ndim; i++) {
      if (y[i] != zero) {
        info = MatSetValues(H,1,&i,1,&j,&y[i],ADD_VALUES); CHKERRQ(info);
      }
    }
    info = VecRestoreArray(user->y,&y); CHKERRQ(info);

  }

  *flg=SAME_NONZERO_PATTERN;

  /* Assemble matrix  */
  info = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);

  return 0;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MatrixFreeHessian"
/* 
   MatrixFreeHessian - Sets a pointer for use in computing Hessian-vector
   products.
    
   Input Parameters:
.  taoapp - the TAO_APPLICATION context
.  X    - the input vector
.  ptr  - optional user-defined context, as set by TaoSetHessian()
   
   Output Parameters:
.  H     - Hessian matrix
.  PrecH - optionally different preconditioning Hessian
.  flag  - flag indicating matrix structure
*/
int MatrixFreeHessian(TAO_APPLICATION taoapp,Vec X,Mat *H,Mat *PrecH,
                      MatStructure *flag,void *ptr)
{
  AppCtx     *user = (AppCtx *) ptr;

  /* Sets location of vector for use in computing matrix-vector products
     of the form H(X)*y  */

  user->xvec = X;   
  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "HessianProductMat"
/* 
   HessianProductMat - Computes the matrix-vector product
   y = mat*svec.

   Input Parameters:
.  mat  - input matrix
.  svec - input vector

   Output Parameters:
.  y    - solution vector
*/
int HessianProductMat(Mat mat,Vec svec,Vec y)
{
  void *ptr;
  MatShellGetContext(mat,&ptr);
  HessianProduct(ptr,svec,y);

  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "HessianProduct"
/* 
   Hessian Product - Computes the matrix-vector product: 
   y = f''(x)*svec.

   Input Parameters
.  ptr  - pointer to the user-defined context
.  svec - input vector

   Output Parameters:
.  y    - product vector
*/
int HessianProduct(void *ptr,Vec svec,Vec y)
{
  AppCtx *user = (AppCtx *)ptr;
  PetscScalar p5 = 0.5, zero = 0.0, one = 1.0, hx, hy, val, area, *x, *s;
  double v, vb, vl, vr, vt, hxhx, hyhy;
  int    info;
  PetscInt nx, ny, i, j, k, ind;

  nx   = user->mx;
  ny   = user->my;
  hx   = user->hx;
  hy   = user->hy;
  hxhx = one/(hx*hx);
  hyhy = one/(hy*hy);

  /* Get pointers to vector data */
  info = VecGetArray(user->xvec,&x); CHKERRQ(info);
  info = VecGetArray(svec,&s); CHKERRQ(info);

  /* Initialize product vector to zero */
  info = VecSet(y, zero); CHKERRQ(info);

  /* Compute f''(x)*s over the lower triangular elements */
  for (j=-1; j<ny; j++) {
    for (i=-1; i<nx; i++) {
       k = nx*j + i;
       v = zero;
       vr = zero;
       vt = zero;
       if (i != -1 && j != -1) v = s[k];
       if (i != nx-1 && j != -1) {
         vr = s[k+1];
         ind = k+1; val = hxhx*(vr-v);
         info = VecSetValues(y,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
       }
       if (i != -1 && j != ny-1) {
         vt = s[k+nx];
         ind = k+nx; val = hyhy*(vt-v);
         info = VecSetValues(y,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
       }
       if (i != -1 && j != -1) {
         ind = k; val = hxhx*(v-vr) + hyhy*(v-vt);
         info = VecSetValues(y,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
       }
     }
   }
  
  /* Compute f''(x)*s over the upper triangular elements */
  for (j=0; j<=ny; j++) {
    for (i=0; i<=nx; i++) {
       k = nx*j + i;
       v = zero;
       vl = zero;
       vb = zero;
       if (i != nx && j != ny) v = s[k];
       if (i != nx && j != 0) {
         vb = s[k-nx];
         ind = k-nx; val = hyhy*(vb-v);
         info = VecSetValues(y,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
       }
       if (i != 0 && j != ny) {
         vl = s[k-1];
         ind = k-1; val = hxhx*(vl-v);
         info = VecSetValues(y,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
       }
       if (i != nx && j != ny) {
         ind = k; val = hxhx*(v-vl) + hyhy*(v-vb);
         info = VecSetValues(y,1,&ind,&val,ADD_VALUES); CHKERRQ(info);
       }
    }
  }
  /* Restore vector data */
  info = VecRestoreArray(svec,&s); CHKERRQ(info);
  info = VecRestoreArray(user->xvec,&x); CHKERRQ(info);

  /* Assemble vector */
  info = VecAssemblyBegin(y); CHKERRQ(info);
  info = VecAssemblyEnd(y); CHKERRQ(info);
 
  /* Scale resulting vector by area */
  area = p5*hx*hy;
  info = VecScale(y, area); CHKERRQ(info);

  info = PetscLogFlops(nx*ny*18); CHKERRQ(info);
  
  return 0;
}


