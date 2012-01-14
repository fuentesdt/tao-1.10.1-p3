/* Chebyshev Quadrature (multiprocessor)
   /home/grignon/lstsqrs/inprogress/chebyq.c
   makefile - /home/grignon/lstsqrs/inprogress/makefile
*/



#include "petscda.h"
#include "tao.h"
#include <math.h>

/* Program usage: chebyq1 [-help] [all TAO options] */


static char help[] = 
"This example demonstrates use of the TAO package to solve a nonlinear \n\
least squares problem.  This example is based on a problem from the\n\
MINPACK-2 test suit.  The Chebyshev problem arises from the determination\n\
of the nodes of a quadature formula with equal weights. \n\n";

/* T
   Concepts: TAO - Solving a bound constrained minimization problem
   Routines: TaoInitialize(); TaoFinalize(); TaoCreate(); TaoSetMethod();
   Routines: TaoSetPetscJacobian(); 
   Routines: TaoSetPetscConstraintsFunction();
   Routines: TaoSetPetscInitialVector(); TaoSetPetscVariableBounds(), TaoSolve();
   Routines: TaoView(); TaoDestroy();
   Processors: N
T */



/* User-defined application context */
typedef struct {
  int   nvar;      /* number of variables */
  int   ncnst;     /* number of constraints */
  int   *row, *column;
  double  *fvec;
  /*Vec   myX;*/
  /*VecScatter scatter;*/
} AppCtx;

/* User provided Routines */
int InitializeVectors(Vec, Vec);
int FormStartingPoint(AppCtx *user, Vec);
int EvaluateConstraints(TAO_SOLVER, Vec, Vec, void *user);
int EvaluateJacobian(TAO_SOLVER, Vec , Mat*, void *ptr);
int AppCtxInitialize(AppCtx*);
int AppCtxDestroy(AppCtx*);

/*--------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **argv)
{

  AppCtx     user;               /* user-defined work context */
  int        info;
  int        m,n;
  double     zero = 0.0, one = 1.0;
  Vec        X, F, G;            /* solution, function, gradient vectors*/
  Vec        XL, XU;             /* lower and upper bound vectors */
  Mat        J;                  /* Jacobian matrix */ 
  IS         is;
  TAO_SOLVER tao;                /* TAO_SOLVER solver context */
  TAO_APPLICATION taoapp;
  TaoTruth success;

   /* Initialize TAO */
  PetscInitialize(&argc,&argv,(char *)0,help);
  TaoInitialize(&argc,&argv,(char *)0,help);


  /* Inititalize user-defined application context */
  info = AppCtxInitialize(&user); CHKERRQ(info);

  /* Allocate vectors */
  info = VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,user.nvar,&X); 
  CHKERRQ(info);
  info = VecDuplicate(X, &XL);CHKERRQ(info);
  info = VecDuplicate(X, &XU);CHKERRQ(info);
  info = VecDuplicate(X, &G);CHKERRQ(info);
  info = VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,user.ncnst,&F); 
  CHKERRQ(info);

  /*
  info = VecCreateSeq(MPI_COMM_SELF,user.nvar,&myX); CHKERRQ(info);
  info = ISCreateStride(MPI_COMM_SELF,user.nvar,0,1,&is); CHKERRQ(info);
  info = VecScatterCreate(X, is, user.myX, &scatter); CHKERRQ(info);
  */

  info = VecGetLocalSize(F,&m);CHKERRQ(info);
  info = VecGetLocalSize(X,&n);CHKERRQ(info);

  /* Create the Jacobian matrix. */
  /*  info =  MatCreateMPIDense(MPI_COMM_WORLD,m,n,user.ncnst,user.nvar,TAO_NULL,&J); */
  info = MatCreateMPIAIJ(MPI_COMM_WORLD,m,n,user.ncnst,user.nvar,user.nvar,PETSC_NULL,user.nvar,PETSC_NULL,&J);  
  CHKERRQ(info);

  /* Set a lower and an upper bound for x. */
  info = VecSet(XL, zero); CHKERRQ(info);
  info = VecSet(XU, one); CHKERRQ(info);

  /* Create TAO solver and set desired solution method */
  info = TaoPetscApplicationCreate(MPI_COMM_WORLD,&taoapp); CHKERRQ(info);

  info = TaoCreate(MPI_COMM_WORLD,"tao_nlsq",&tao);CHKERRQ(info);

  /* Set the function and Jacobian routines. */
  info = TaoSetPetscFunctionGradient(taoapp,X,G,TAO_NULL,TAO_NULL); CHKERRQ(info);
  info = TaoSetPetscJacobian(taoapp, J, EvaluateJacobian, &user);  CHKERRQ(info);
  info = TaoSetPetscConstraintsFunction(taoapp, F, EvaluateConstraints, &user); 
  CHKERRQ(info);

  /*  Compute the standard starting point. */
  info = FormStartingPoint(&user, X); CHKERRQ(info);
  info = TaoSetPetscInitialVector(taoapp,G); CHKERRQ(info);
  info = TaoSetPetscVariableBounds(taoapp,XL,XU);CHKERRQ(info);

  info = TaoSetApplication(tao,taoapp); CHKERRQ(info); 
  info = TaoSetFromOptions(tao); CHKERRQ(info);

  /* Perform the Solve */
  info = TaoSolve(tao,&success); CHKERRQ(info);


  /* View the vectors and matrices to check for correctness */
  /* info = VecView(F,VIEWER_STDOUT_SELF); CHKERRQ(info);
     info = MatView(J,VIEWER_STDOUT_SELF); CHKERRQ(info); */

  /* View TAO solver information */
  info = TaoView(tao); CHKERRQ(info);


   /* Free data structures */
  info = VecDestroy(X); CHKERRQ(info);
  info = VecDestroy(G); CHKERRQ(info);
  info = VecDestroy(XL); CHKERRQ(info);
  info = VecDestroy(XU); CHKERRQ(info);
  info = VecDestroy(F); CHKERRQ(info);
  info = MatDestroy(J); CHKERRQ(info);
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoApplicationDestroy(taoapp); CHKERRQ(info);
  /*  info = ISDestroy(is); CHKERRQ(info); */

  /* Destroy user-defined application context */
  info = AppCtxDestroy(&user); CHKERRQ(info);

  /* Finalize TAO */
  TaoFinalize();
  PetscFinalize();
  return 0;     
}



/*--------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "FormStartingPoint"
int FormStartingPoint(AppCtx *ptr, Vec X){

  AppCtx *user = (AppCtx *)ptr;
  int        i, low,high,info;
  double     dx;
  double     *x;         
  double     one = 1.0;

  dx = one/(user->nvar+one);
  info = VecGetOwnershipRange(X,&low,&high); CHKERRQ(info);

  info = VecGetArray(X,&x); CHKERRQ(info);
  for(i = low; i < high; i++){
    x[i] = i*dx;
  }
  info = VecRestoreArray(X,&x); CHKERRQ(info);

  return 0;
}


/*--------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "EvaluateConstraints"
int EvaluateConstraints(TAO_SOLVER tao, Vec X, Vec F, void *ptr){

  AppCtx *user = (AppCtx *)ptr;
  int        i,j, iev, info;
  int        low, high;
  double     zero=0.0, one=1.0, two=2.0;
  double     ti, temp, temp1, temp2;
  double     *xvec;
  double     dx;
  
  dx = one/(user->nvar);

  /*
  info = VecScatterBegin(X,user->myX,INSERT_VALUES, SCATTER_FORWARD, 
  user->scatter); CHKERRQ(info);
  info = VecScatterEnd(X,user->myX,INSERT_VALUES, SCATTER_FORWARD, 
  user->scatter); CHKERRQ(info);
  */

  info = VecGetArray(X,&xvec); CHKERRQ(info);
  info = VecGetOwnershipRange(X,&low,&high);CHKERRQ(info);
  info = VecSet(F, zero); CHKERRQ(info);

  for(i=0; i<user->ncnst; i++){
    user->fvec[i] = zero;
  }

  for(j=low; j<high; j++){
    temp1 = one;
    temp2 = two*xvec[j-low] - one;
    temp = two*temp2;
    for(i=0;i<user->ncnst;i++){
      user->fvec[i] = user->fvec[i] + temp2;
      ti = temp*temp2 - temp1;
      temp1 = temp2;
      temp2 = ti;
    }
  }

  iev = -1;
  for(i=0; i<user->ncnst; i++){
    user->fvec[i-low] = dx*user->fvec[i];
    if (iev > 0) user->fvec[i] = user->fvec[i] + one/(((i+1)*(i+1))-one);
    iev = -iev;
  }


  info = VecRestoreArray(X,&xvec); CHKERRQ(info);
  info = VecSetValues(F,user->ncnst,user->row,user->fvec,ADD_VALUES); 
  CHKERRQ(info);

  info = VecAssemblyBegin(F);
  info = VecAssemblyEnd(F);
  

  return 0;

}


/*--------------------------------------------------------------------*/
#undef __FUNC__
#define __FUNC__ "EvaluateJacobian"
int EvaluateJacobian(TAO_SOLVER tao, Vec X, Mat *JJ, void *ptr){

  AppCtx     *user = (AppCtx *)ptr;
  int        i,j, info;
  int        low, high;
  double     zero=0.0, one=1.0, two=2.0, four=4.0;
  double     ti, temp1, temp2, temp3, temp4;
  double     twoval,val;
  double     *xvec;
  double     dx;
  double     *jaccol;
  Mat J=*JJ;

  info =PetscMalloc(sizeof(double)*user->ncnst,&jaccol); CHKERRQ(info);
  CHKERRQ(info); 

  dx = one/user->nvar;
  info = VecGetOwnershipRange(X,&low,&high);CHKERRQ(info);
  info = VecGetArray(X,&xvec); CHKERRQ(info);

  for(j=low; j<high; j++){
    val =  two*xvec[j-low] - one;
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

  info = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);


  info = PetscFree(jaccol); CHKERRQ(info);


  /*  MatView(*J, VIEWER_STDOUT_WORLD); */
  
  info = VecRestoreArray(X,&xvec); CHKERRQ(info);
  return 0;

}




/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "AppCtxInitialize"
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

  PetscFree(user->row);
  PetscFree(user->column);
  PetscFree(user->fvec);

  return 0;

}

