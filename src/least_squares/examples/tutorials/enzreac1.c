/*$Id$*/

/* Program usage: enzreac [-help] [all TAO options] */



/* 
   Include "tao.h" so that we can use TAO solvers.  Note that this
   file automatically includes libraries such as:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - sysem routines        petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

*/

#include "tao.h"

static char help[] = 
"Demonstrates the use of the TAO package to solve\n\
nonlinear least squares problme.  This example is based on a problem\n\
from the MINPACK-2 test suite. This problem arises from the analysis of\n\
the kinetic data for an enzyme reaction. \n\n";

/* T
   Concepts: TAO - Solving a system of nonlinear equations, nonlinear least squares;
   Routines: TaoInitialize(); TaoFinalize(); 
   Routines: TaoCreate();  TaoDestroy();
   Routines: TaoPetscApplicationCreate(); TaoApplicationDestroy();
   Routines: TaoSetPetscFunction(); TaoSetPetscJacobian(); 
   Routines: TaoSetPetscConstraintsFunction();
   Routines: TaoSetPetscInitialVector(); TaoSetApplication();
   Routines: TaoSolve(); TaoSetFromOptions(); 
   Processors: 1
T*/ 

const int nvars = 4;
const int nconsts = 11;

/* User-defined application context */
typedef struct {
  Vec V, Y;
} AppCtx;


int InitializeVectors(Vec, Vec);
int FormStartingPoint(Vec);
int EvaluateConstraints(TAO_SOLVER, Vec, Vec, void *);
int EvaluateJacobian(TAO_SOLVER, Vec, Mat*, void *);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  int        info;             /* Used to check for nonzero function returns */
  Vec        x,f;            /* solution, function vectors*/
  Mat        J;                /* Jacobian matrix */
  TAO_SOLVER tao;              /* TAO_SOLVER solver context */
  TAO_APPLICATION taoapp;      /* PETSc application */
  AppCtx     user;             /* application context */
  

  /* Initialize TAO and PETSc */
  PetscInitialize(&argc,&argv,(char *)0,help);
  TaoInitialize(&argc,&argv,(char *)0,help);

  /* Allocate vectors */
  info = VecCreateSeq(MPI_COMM_SELF,nvars,&x); CHKERRQ(info);
  info = VecCreateSeq(MPI_COMM_SELF,nconsts,&f); CHKERRQ(info);
  info = VecDuplicate(f,&user.Y); CHKERRQ(info);
  info = VecDuplicate(f,&user.V); CHKERRQ(info);
  
  /* Initialize application data -- user function */
  info = InitializeVectors(user.V, user.Y); CHKERRQ(info);

  /* Create Jacobian matrix  */
  info = MatCreateSeqDense(MPI_COMM_SELF,nconsts,nvars,TAO_NULL,&J); CHKERRQ(info);

  /* TAO code begines here */

  /* Create TAO solver and set desired solution method */
  info = TaoCreate(MPI_COMM_SELF,"tao_nlsq",&tao);CHKERRQ(info);
  info = TaoPetscApplicationCreate(MPI_COMM_SELF,&taoapp); CHKERRQ(info);


  /* Set the function and Jacobian routines. */
  info = TaoSetPetscFunction(taoapp,x,TAO_NULL,TAO_NULL);  CHKERRQ(info);
  info = TaoSetPetscJacobian(taoapp, J, EvaluateJacobian, (void*)&user); CHKERRQ(info);
  info = TaoSetPetscConstraintsFunction(taoapp, f, EvaluateConstraints, (void*)&user);  CHKERRQ(info);

  /*  Compute the starting point. */
  info = FormStartingPoint(x); CHKERRQ(info);
  info = TaoSetPetscInitialVector(taoapp,x); CHKERRQ(info); 

  /* Now that the PETSc application is set, attach to TAO context */
  info = TaoSetApplication(tao,taoapp); CHKERRQ(info); 

  /* Check for any TAO command line arguments */
  info = TaoSetFromOptions(tao); CHKERRQ(info);

  /* Perform the Solve */
  info = TaoSolve(tao); CHKERRQ(info);

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoApplicationDestroy(taoapp); CHKERRQ(info);

  /* Free PETSc data structures */
  info = VecDestroy(user.V); CHKERRQ(info);
  info = VecDestroy(user.Y); CHKERRQ(info);
  info = VecDestroy(x); CHKERRQ(info);
  info = VecDestroy(f); CHKERRQ(info);
  info = MatDestroy(J); CHKERRQ(info);
   
  /* Finalize TAO */
  TaoFinalize();
  PetscFinalize();

  return 0;     
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateConstraints"
int EvaluateConstraints(TAO_SOLVER tao, Vec X, Vec F, void *ptr){

  AppCtx *user = (AppCtx *)ptr;
  int i, info;
  double temp1, temp2, val;
  double *v, *y, *x;
  int dimx,dimv,dimy;

  info = VecGetArray(X,&x);CHKERRQ(info);
  info = VecGetArray(user->V,&v);CHKERRQ(info);
  info = VecGetArray(user->Y,&y);CHKERRQ(info);
  info = VecGetLocalSize(user->V,&dimv); CHKERRQ(info);
  info = VecGetLocalSize(user->Y,&dimy); CHKERRQ(info);
  info = VecGetLocalSize(X,&dimx); CHKERRQ(info);

  for (i = 0; i < nconsts; i ++){

    temp1 = v[i]*(v[i]+x[1]);
    temp2 = v[i]*(v[i]+x[2]) + x[3];

    val = y[i] - x[0]*temp1/temp2;	 
    info = VecSetValues(F,1,&i,&val,INSERT_VALUES); CHKERRQ(info);
  }

  info = VecRestoreArray(X,&x);CHKERRQ(info);
  info = VecRestoreArray(user->V,&v);CHKERRQ(info);
  info = VecRestoreArray(user->Y,&y);CHKERRQ(info);

  /* Assemble function vector */
  info = VecAssemblyBegin(F); CHKERRQ(info);
  info = VecAssemblyEnd(F); CHKERRQ(info);

  PetscLogFlops(nconsts*8);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateJacobian"
int EvaluateJacobian(TAO_SOLVER tao, Vec X, Mat *JJ, void *ptr){

  AppCtx *user = (AppCtx *) ptr;
  int i, info;
  double temp1, temp2;
  double *v, *x;
  double val[nconsts][nvars];
  int    row[nconsts], col[nvars];
  Mat J=*JJ;

  info = VecGetArray(user->V,&v); CHKERRQ(info);
  info = VecGetArray(X,&x); CHKERRQ(info);

  for (i = 0; i < nconsts; i ++){

    temp1 = v[i]*(v[i]+x[1]);
    temp2 = v[i]*(v[i]+x[2]) + x[3];
    row[i] = i;
    if (i < nvars) col[i] = i;
    
    val[i][0] = -temp1/temp2;
    val[i][1] = -v[i]*x[0]/temp2;
    val[i][2] = val[i][0]*val[i][1];
    val[i][3] = val[i][2]/v[i];
  }

  info = MatSetValues(J,nconsts,row,nvars,col,*val,INSERT_VALUES);  CHKERRQ(info);

  info = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(info)
  info = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);


  info = VecRestoreArray(user->V,&v); CHKERRQ(info);
  info = VecRestoreArray(X,&x); CHKERRQ(info);

  PetscLogFlops(nconsts*10);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "FormStartingPoint"
int FormStartingPoint(Vec X)
{
  double *x;
  int info;

  info = VecGetArray(X,&x); CHKERRQ(info);

  x[0] = 2.5e-1;
  x[1] = 3.9e-1;
  x[2] = 4.15e-1;
  x[3] = 3.9e-1;

  info = VecRestoreArray(X,&x); CHKERRQ(info);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "InitializeVectors"
int InitializeVectors(Vec V, Vec Y) 
{
  double *v;
  double *y;
  int info;

  info = VecGetArray(V, &v); CHKERRQ(info);

  v[0] = 4.0e0;
  v[1] = 2.0e0;
  v[2] = 1.0e0;
  v[3] = 5.0e-1;
  v[4] = 2.5e-1;
  v[5] = 1.67e-1;
  v[6] = 1.25e-1;
  v[7] = 1.0e-1;
  v[8] = 8.33e-2; 
  v[9] = 7.14e-2;
  v[10] = 6.25e-2;

  /* Restore vector data */
  info = VecRestoreArray(V, &v); CHKERRQ(info);

  info = VecGetArray(Y, &y); CHKERRQ(info);

  y[0] = 1.957e-1;
  y[1] = 1.947e-1;
  y[2] = 1.735e-1;
  y[3] = 1.6e-1;
  y[4] = 8.44e-2;
  y[5] = 6.27e-2;
  y[6] = 4.56e-2;
  y[7] = 3.42e-2;
  y[8] = 3.23e-2;
  y[9] = 2.35e-2;
  y[10] = 2.46e-2;

  /* Restore vector data */
  info = VecRestoreArray(Y, &y); CHKERRQ(info);

  return 0;
}

