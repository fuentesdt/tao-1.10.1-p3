/* Program usage: chebyq1 [-help] [all TAO options] */

/* 
   Include "petscda.h" so that we can use distributed arrays (DAs).
   Include "tao.h" so that we can use TAO solvers.  Note that this
   file automatically includes libraries such as:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - sysem routines        petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

*/

#include "petscda.h"
#include "tao.h"


static char help[] = 
"This example demonstrates use of the TAO package to solve a nonlinear \n\
least squares problem.  This example is based on a problem from the\n\
MINPACK-2 test suit.  The Chebyshev problem arises from the determination\n\
of the nodes of a quadature formula with equal weights. \n\n";

/* T
   Concepts: TAO - Solving a system of nonlinear equations, nonlinear least squares;
   Routines: TaoInitialize(); TaoFinalize(); 
   Routines: TaoCreate(); TaoDestroy();
   Routines: TaoPetscApplicationCreate(); TaoApplicationDestroy();
   Routines: TaoSetPetscFunctionGradient(); 
   Routines: TaoSetPetscConstraintsFunction(); TaoSetPetscJacobian(); 
   Routines: TaoSetPetscInitialVector(); TaoSetPetscVariableBounds(); 
   Routines: TaoSetApplication(); TaoSetFromOptions(); TaoSolve(); TaoView(); 
   Processors: 1
T*/



/* User-defined application context */
typedef struct {
  /* Application parameters */
  int   nvar;      /* number of variables */
  int   ncnst;     /* number of constraints */
  
  /* working space */
  int   *row, *column;
  double  *fvec;
} AppCtx;

/* User provided Routines */
int InitializeVectors(Vec, Vec);
int FormStartingPoint(AppCtx *user, Vec);
int EvaluateConstraints(TAO_SOLVER, Vec, Vec, void *);
int EvaluateJacobian(TAO_SOLVER, Vec, Mat*, void *);
int AppCtxInitialize(AppCtx*);
int AppCtxDestroy(AppCtx*);


/*--------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  int        info;               /* used to check for functions returning nonzeros */
  Vec        x, f;               /* solution, function */
  Vec        xl, xu;             /* lower and upper bound vectors */
  Mat        J;                  /* Jacobian matrix */ 
  TAO_SOLVER tao;                /* TAO_SOLVER solver context */
  TAO_APPLICATION taoapp;        /* The PETSc application */
  int        iter;               /* iteration information */
  double     ff;
  double     zero = 0.0, one = 1.0;
  AppCtx     user;               /* user-defined work context */

   /* Initialize TAO and PETSc */
  PetscInitialize(&argc,&argv,(char *)0,help);
  TaoInitialize(&argc,&argv,(char *)0,help);


  /* Inititalize user-defined application context */
  info = AppCtxInitialize(&user); CHKERRQ(info);

  /* Allocate vectors */
  info = VecCreateSeq(MPI_COMM_SELF,user.nvar,&x); CHKERRQ(info);
  info = VecCreateSeq(MPI_COMM_SELF,user.ncnst,&f); CHKERRQ(info);


  /* Create the Jacobian matrix. */
  info = MatCreateSeqDense(MPI_COMM_SELF,user.ncnst,user.nvar,PETSC_NULL,&J);  
      CHKERRQ(info);


  /* TAO code begins here */

  /* Create TAO solver and set desired solution method */
  info = TaoCreate(MPI_COMM_SELF,"tao_nlsq",&tao);CHKERRQ(info);
  info = TaoPetscApplicationCreate(MPI_COMM_SELF,&taoapp); CHKERRQ(info);


  /* Set the function and Jacobian routines. */
  info = TaoSetPetscFunction(taoapp,x,TAO_NULL,TAO_NULL); CHKERRQ(info);
  info = TaoSetPetscJacobian(taoapp, J, EvaluateJacobian, (void*)&user);  CHKERRQ(info);
  info = TaoSetPetscConstraintsFunction(taoapp, f, EvaluateConstraints, (void*)&user); CHKERRQ(info);


  /* Set a lower and an upper bound for x. */
  info = VecDuplicate(x, &xl);CHKERRQ(info);
  info = VecDuplicate(x, &xu);CHKERRQ(info);
  info = VecSet(xl, zero); CHKERRQ(info);
  info = VecSet(xu, one); CHKERRQ(info);
  info = TaoSetPetscVariableBounds(taoapp,xl,xu);CHKERRQ(info);

  /*  Compute the starting point. */
  info = FormStartingPoint(&user, x); CHKERRQ(info);
  info = TaoSetPetscInitialVector(taoapp,x); CHKERRQ(info);

  /* Now that the PETSc application is set, attach to TAO context */
  info = TaoSetApplication(tao,taoapp); CHKERRQ(info); 

  /* Check for any TAO command line arguments */
  info = TaoSetFromOptions(tao); CHKERRQ(info);

  /* Perform the Solve */
  info = TaoSolve(tao); CHKERRQ(info);

  /* View iteration data */
  info = TaoGetIterationData(tao,&iter,&ff,0,0,0,0); CHKERRQ(info);
  PetscPrintf(PETSC_COMM_SELF,"Solved: Iterations: %d, Residual: %5.3e\n",
	      iter,ff);

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoApplicationDestroy(taoapp); CHKERRQ(info);

   /* Free PETSc data structures */
  info = VecDestroy(x); CHKERRQ(info);
  info = VecDestroy(xl); CHKERRQ(info);
  info = VecDestroy(xu); CHKERRQ(info);
  info = VecDestroy(f); CHKERRQ(info);
  info = MatDestroy(J); CHKERRQ(info);

  /* Destroy user-defined application context */
  info = AppCtxDestroy(&user); CHKERRQ(info);

  /* Finalize TAO */
  TaoFinalize();
  PetscFinalize();
  return 0;     
}



/*--------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "FormStartingPoint"
int FormStartingPoint(AppCtx *ptr, Vec X){

  AppCtx *user = (AppCtx *)ptr;
  int        i, info;
  double     dx;
  double     *x;         
  double     one = 1.0;

  dx = one/(user->nvar+one);

  info = VecGetArray(X,&x); CHKERRQ(info);
  for(i = 0; i < user->nvar; i++){
    x[i] = i*dx;
  }
  info = VecRestoreArray(X,&x); CHKERRQ(info);

  return 0;
}


/*--------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "EvaluateConstraints"
int EvaluateConstraints(TAO_SOLVER tao, Vec X, Vec Res, void *ptr){

  AppCtx *user = (AppCtx *)ptr;
  int        i,j, iev, info;
  double     zero=0.0, one=1.0, two=2.0;
  double     ti, temp, temp1, temp2;
  double     *xvec,*fvec;
  double     dx;

  dx = one/(user->nvar);
  info = VecGetArray(X,&xvec); CHKERRQ(info);
  info = VecGetArray(Res,&fvec); CHKERRQ(info);

  for(i=0; i<user->ncnst; i++){
    fvec[i] = zero;
  }

  for(j=0; j<user->nvar; j++){
    temp1 = one;
    temp2 = two*xvec[j] - one;
    temp = two*temp2;
    for(i=0;i<user->ncnst;i++){
      fvec[i] = fvec[i] + temp2;
      ti = temp*temp2 - temp1;
      temp1 = temp2;
      temp2 = ti;
    }
  }

  iev = -1;
  for(i=0; i<user->ncnst; i++){
    fvec[i] = dx*fvec[i];
    if (iev > 0) fvec[i] = fvec[i] + one/(((i+1)*(i+1))-one);
    iev = -iev;
  }

  PetscLogFlops(user->nvar * (3 + user->ncnst*3) + user->ncnst);
  info = VecRestoreArray(X,&xvec); CHKERRQ(info);
  info = VecRestoreArray(Res,&fvec); CHKERRQ(info);

  return 0;

}


/*--------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "EvaluateJacobian"
int EvaluateJacobian(TAO_SOLVER tao, Vec X, Mat *JJ, void *ptr){

  AppCtx     *user = (AppCtx *)ptr;
  int        i,j, info;
  double     zero=0.0, one=1.0, two=2.0, four=4.0;
  double     ti, temp1, temp2, temp3, temp4;
  double     twoval,val;
  double     *xvec;
  double     dx;
  double     *jaccol=user->fvec;
  Mat J=*JJ;

  dx = one/user->nvar;
  info = VecGetArray(X,&xvec); CHKERRQ(info);

  for(j=0; j<user->nvar; j++){
    val =  two*xvec[j] - one;
    twoval = two*val;

    temp1 = one;
    temp2 = val;
    temp3 = zero;
    temp4 = two;

    for(i=0; i<user->ncnst; i++){
      jaccol[i] = dx*temp4;
      ti = four*temp2 + twoval*temp4 - temp3;
      temp3 = temp4;
      temp4 = ti;
      ti = twoval*temp2 - temp1;
      temp1 = temp2;
      temp2 = ti;
      
    }
    info = MatSetValues(J,user->ncnst, user->row, 1,&j,jaccol,INSERT_VALUES);  CHKERRQ(info);   

  }

  PetscLogFlops(user->nvar * (3 + 7*user->ncnst));
  info = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);

  info = VecRestoreArray(X,&xvec); CHKERRQ(info);
  return 0;

}




/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "AppCtxInitialize"
/*
    AppCtxInitialize - Sets up user-defined work space; initializes problem
    parameters.

    Input Parameter:
    user - user-defined application context

    Output Parameter:
    user - user-defined application context
*/
int AppCtxInitialize(AppCtx *user)
{
  int      j, info;
  PetscTruth flg;

  user->ncnst = 8; user->nvar = 8;
  info = PetscOptionsGetInt(TAO_NULL,"-ncnst",&user->ncnst,&flg); 
  CHKERRQ(info);
  info = PetscOptionsGetInt(TAO_NULL,"-nvar",&user->nvar,&flg); 
  CHKERRQ(info);

  info = PetscMalloc(sizeof(int)*user->ncnst,&user->row);CHKERRQ(info);
  for(j=0; j < user->ncnst; j++){
    user->row[j] = j;
  }

  info = PetscMalloc(sizeof(int)*user->nvar,&user->column); CHKERRQ(info);
  for(j=0; j < user->nvar; j++){
    user->column[j] = j;
  }

  info = PetscMalloc(sizeof(double)*user->ncnst,&user->fvec); CHKERRQ(info);

  return 0;
}

int AppCtxDestroy(AppCtx *user){
  int info;

  info = PetscFree(user->row); CHKERRQ(info);
  info = PetscFree(user->column); CHKERRQ(info);
  info = PetscFree(user->fvec); CHKERRQ(info);

  return 0;

}

