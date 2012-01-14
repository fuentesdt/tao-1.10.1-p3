/**********************************************************************
    American Put Options Pricing using the Black-Scholes Equation
  
   Background (European Options):
     The standard European option is a contract where the holder has the right
     to either buy (call option) or sell (put option) an underlying asset at 
     a designated future time and price.
  
     The classic Black-Scholes model begins with an assumption that the
     price of the underlying asset behaves as a lognormal random walk.  
     Using this assumption and a no-arbitrage argument, the following
     linear parabolic partial differential equation for the value of the
     option results:
  
       dV/dt + 0.5(sigma**2)(S**alpha)(d2V/dS2) + (r - D)S(dV/dS) - rV = 0.
  
     Here, sigma is the volatility of the underling asset, alpha is a
     measure of elasticity (typically two), D measures the dividend payments
     on the underling asset, and r is the interest rate.
  
     To completely specify the problem, we need to impose some boundary
     conditions.  These are as follows:
  
       V(S, T) = max(E - S, 0)
       V(0, t) = E for all 0 <= t <= T
       V(s, t) = 0 for all 0 <= t <= T and s->infinity
  
     where T is the exercise time time and E the strike price (price paid
     for the contract).
  
     An explicit formula for the value of an European option can be
     found.  See the references for examples.
  
   Background (American Options):
     The American option is similar to its European counterpart.  The
     difference is that the holder of the American option can excercise
     their right to buy or sell the asset at any time prior to the
     expiration.  This additional ability introduce a free boundary into
     the Black-Scholes equation which can be modeled as a linear
     complementarity problem.
  
       0 <= -(dV/dt + 0.5(sigma**2)(S**alpha)(d2V/dS2) + (r - D)S(dV/dS) - rV)
         complements 
       V(S,T) >= max(E-S,0)
  
     where the variables are the same as before and we have the same boundary
     conditions.
  
     There is not explicit formula for calculating the value of an American
     option.  Therefore, we discretize the above problem and solve the
     resulting linear complementarity problem.
  
     We will use backward differences for the time variables and central
     differences for the space variables.  Crank-Nicholson averaging will
     also be used in the discretization.  The algorithm used by the code
     solves for V(S,t) for a fixed t and then uses this value in the
     calculation of V(S,t - dt).  The method stops when V(S,0) has been
     found.
  
   References:
     Huang and Pang, "Options Pricing and Linear Complementarity,"
       Journal of Computational Finance, volume 2, number 3, 1998.
     Wilmott, "Derivatives: The Theory and Practice of Financial Engineering,"
       John Wiley and Sons, New York, 1998.
***************************************************************************/

/*
  Include "tao.h" so we can use TAO solvers with PETSc support.  
  Include "petscda.h" so that we can use distributed arrays (DAs) for managing
  the parallel mesh.
*/

#include "petscda.h"
#include "tao.h"

static char  help[] =
"This example demonstrates use of the TAO package to\n\
solve a linear complementarity problem for pricing American put options.\n\
The code uses backward differences in time and central differences in\n\
space.  The command line options are:\n\
  -rate <r>, where <r> = interest rate\n\
  -sigma <s>, where <s> = volatility of the underlying\n\
  -alpha <a>, where <a> = elasticity of the underlying\n\
  -delta <d>, where <d> = dividend rate\n\
  -strike <e>, where <e> = strike price\n\
  -expiry <t>, where <t> = the expiration date\n\
  -mt <tg>, where <tg> = number of grid points in time\n\
  -ms <sg>, where <sg> = number of grid points in space\n\
  -es <se>, where <se> = ending point of the space discretization\n\n";

/*T
   Concepts: TAO - Solving a complementarity problem
   Routines: TaoInitialize(); TaoFinalize();
   Routines: TaoCreate(); TaoDestroy();
   Routines: TaoApplicationCreate(); TaoAppDestroy();
   Routines: TaoAppSetJacobianRoutine(); TaoAppSetConstraintRoutine();
   Routines: TaoAppSetJacobianMat(); TaoSetOptions();
   Routines: TaoSolveApplication();
   Routines: TaoAppSetVariableBoundsRoutine(); TaoAppSetInitialSolutionVec();
   Processors: 1
T*/


/*
  User-defined application context - contains data needed by the
  application-provided call-back routines, FormFunction(), and FormJacobian().
*/

typedef struct {
  double *Vt1;                /* Value of the option at time T + dt */
  double *c;                  /* Constant -- (r - D)S */
  double *d;	              /* Constant -- -0.5(sigma**2)(S**alpha) */

  PetscReal rate;                /* Interest rate */
  PetscReal sigma, alpha, delta; /* Underlying asset properties */
  PetscReal strike, expiry;      /* Option contract properties */

  PetscReal es;		      /* Finite value used for maximum asset value */
  PetscReal ds, dt;	      /* Discretization properties */
  PetscInt ms, mt;		      /* Number of elements */

  DA da;
} AppCtx;

/* -------- User-defined Routines --------- */

int FormConstraints(TAO_APPLICATION, Vec, Vec, void *);
int FormJacobian(TAO_APPLICATION, Vec, Mat *, Mat*, MatStructure*, void *);
int ComputeVariableBounds(TAO_APPLICATION, Vec, Vec, void*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  int      info;          /* used to check for functions returning nonzeros */
  Vec      x;             /* solution vector */
  Vec      c;             /* Constraints function vector */
  Mat J;		  /* jacobian matrix */
  PetscTruth flg;	  /* A return variable when checking for user options */
  TAO_SOLVER tao;	  /* TAO_SOLVER solver context */
  TAO_APPLICATION my_app; /* The PETSc application */
  AppCtx user;		  /* user-defined work context */
  PetscInt i, j;               
  PetscInt    xs,xm,gxs,gxm;
  double sval = 0;
  PetscScalar *x_array;
  Vec    localX;

  /* Initialize PETSc, TAO */
  PetscInitialize(&argc, &argv, (char *)0, help);
  TaoInitialize(&argc, &argv, (char *)0, help);

  /* 
     Initialize the user-defined application context with reasonable 
     values for the American option to price 
  */
  user.rate = 0.04;
  user.sigma = 0.40;
  user.alpha = 2.00;
  user.delta = 0.01;
  user.strike = 10.0;
  user.expiry = 1.0;
  user.mt = 10;
  user.ms = 150;
  user.es = 100.0;
  
  /* Read in alternative values for the American option to price */
  info = PetscOptionsGetReal(PETSC_NULL, "-alpha", &user.alpha, &flg); 
         CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL, "-delta", &user.delta, &flg); 
         CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL, "-es", &user.es, &flg); 
         CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL, "-expiry", &user.expiry, &flg); 
         CHKERRQ(info);
  info = PetscOptionsGetInt(PETSC_NULL, "-ms", &user.ms, &flg); 
         CHKERRQ(info);
  info = PetscOptionsGetInt(PETSC_NULL, "-mt", &user.mt, &flg); 
         CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL, "-rate", &user.rate, &flg); 
         CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL, "-sigma", &user.sigma, &flg); 
         CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL, "-strike", &user.strike, &flg); 
         CHKERRQ(info);

  /* Check that the options set are allowable (needs to be done) */

  user.ms++;
  info = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,user.ms,1,1,
                    PETSC_NULL,&user.da); CHKERRQ(info);
  /* Create appropriate vectors and matrices */

  info = DAGetCorners(user.da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(user.da,&gxs,PETSC_NULL,PETSC_NULL,&gxm,PETSC_NULL,PETSC_NULL); CHKERRQ(info);

  info = DACreateGlobalVector(user.da,&x);CHKERRQ(info);
  /* 
     Finish filling in the user-defined context with the values for
     dS, dt, and allocating space for the constants 
  */
  user.ds = user.es / (user.ms-1);
  user.dt = user.expiry / user.mt;

  info = PetscMalloc((gxm)*sizeof(double),&(user.Vt1)); CHKERRQ(info);
  info = PetscMalloc((gxm)*sizeof(double),&(user.c)); CHKERRQ(info);
  info = PetscMalloc((gxm)*sizeof(double),&(user.d)); CHKERRQ(info);

  /* 
     Calculate the values for the constant.  Vt1 begins with the ending 
     boundary condition.  
  */
  for (i=0; i<gxm; i++) {
    sval = (gxs+i)*user.ds;
    user.Vt1[i] = PetscMax(user.strike - sval, 0);
    user.c[i] = (user.delta - user.rate)*sval;
    user.d[i] = -0.5*user.sigma*user.sigma*pow(sval, user.alpha);
  }
  if (gxs+gxm==user.ms){
    user.Vt1[gxm-1] = 0;
  }

  info = VecDuplicate(x, &c); CHKERRQ(info);

  /* 
     Allocate the matrix used by TAO for the Jacobian.  Each row of
     the Jacobian matrix will have at most three elements.
  */
  info = DAGetMatrix(user.da,MATAIJ,&J);CHKERRQ(info);
	 
  /* The TAO code begins here */

  /* Create TAO solver and set desired solution method  */
  info = TaoCreate(PETSC_COMM_WORLD, "tao_ssils", &tao); CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_WORLD, &my_app); CHKERRQ(info);

  /* Set routines for constraints function and Jacobian evaluation */
  info = TaoAppSetConstraintRoutine(my_app, FormConstraints, (void *)&user); 
  CHKERRQ(info);
  info = TaoAppSetJacobianRoutine(my_app, FormJacobian, (void *)&user); CHKERRQ(info); 
  info = TaoAppSetJacobianMat(my_app, J, J); CHKERRQ(info);
    
  /* Set the variable bounds */
  info = TaoAppSetVariableBoundsRoutine(my_app,ComputeVariableBounds,(void*)&user); 
  CHKERRQ(info);

  /* Set initial solution guess */
  info = VecGetArray(x,&x_array); CHKERRQ(info);
  for (i=0; i< xm; i++) 
    x_array[i] = user.Vt1[i-gxs+xs];
  info = VecRestoreArray(x,&x_array); CHKERRQ(info);
  /* Set data structure */ 
  info = TaoAppSetInitialSolutionVec(my_app, x); CHKERRQ(info);

  /* Set routines for function and Jacobian evaluation */
  info = TaoSetOptions(my_app,tao); CHKERRQ(info);

  /* Iteratively solve the linear complementarity problems  */
  for (i = 1; i < user.mt; i++) {

    /* Solve the current version */
    info = TaoSolveApplication(my_app,tao);  CHKERRQ(info);

    /* Update Vt1 with the solution */
    info = DAGetLocalVector(user.da,&localX);CHKERRQ(info);
    info = DAGlobalToLocalBegin(user.da,x,INSERT_VALUES,localX); CHKERRQ(info);
    info = DAGlobalToLocalEnd(user.da,x,INSERT_VALUES,localX); CHKERRQ(info);
    info = VecGetArray(localX,&x_array); CHKERRQ(info);
    for (j = 0; j < gxm; j++) {
      user.Vt1[j] = x_array[j];
    }
    info = VecRestoreArray(x,&x_array); CHKERRQ(info);
    info = DARestoreLocalVector(user.da,&localX); CHKERRQ(info);
  }

  /* Free TAO data structures */

  info = TaoDestroy(tao); CHKERRQ(info);

  info = TaoAppDestroy(my_app);  CHKERRQ(info);

  /* Free PETSc data structures */
  info = VecDestroy(x); CHKERRQ(info);
  info = VecDestroy(c); CHKERRQ(info);
  info = MatDestroy(J); CHKERRQ(info);
  info = DADestroy(user.da); CHKERRQ(info);
  /* Free user-defined workspace */
  info = PetscFree(user.Vt1); CHKERRQ(info);
  info = PetscFree(user.c); CHKERRQ(info);
  info = PetscFree(user.d); CHKERRQ(info);

  /* Finalize TAO and PETSc */
  PetscFinalize();
  TaoFinalize();

  return 0;
}

/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "ComputeVariableBounds"
int ComputeVariableBounds(TAO_APPLICATION tao, Vec xl, Vec xu, void*ctx)
{
  AppCtx *user = (AppCtx *) ctx;
  int info;
  PetscInt  i;
  PetscInt  xs,xm;
  PetscInt  ms = user->ms;
  PetscScalar sval=0.0,*xl_array,ub= TAO_INFINITY;

  /* Set the variable bounds */
  info = VecSet(xu, ub); CHKERRQ(info);
  info = DAGetCorners(user->da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL); CHKERRQ(info);

  info = VecGetArray(xl,&xl_array); CHKERRQ(info);
  for (i = 0; i < xm; i++){
    sval = (xs+i)*user->ds;
    xl_array[i] = PetscMax(user->strike - sval, 0);
  }
  info = VecRestoreArray(xl,&xl_array); CHKERRQ(info);

  if (xs==0){
    info = VecGetArray(xu,&xl_array); CHKERRQ(info);
    xl_array[0] = PetscMax(user->strike, 0);
    info = VecRestoreArray(xu,&xl_array); CHKERRQ(info);
  }
  if (xs+xm==ms){
    info = VecGetArray(xu,&xl_array); CHKERRQ(info);
    xl_array[xm-1] = 0;
    info = VecRestoreArray(xu,&xl_array); CHKERRQ(info);
  }

  return 0;
}
/* -------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormConstraints"

/*  
    FormFunction - Evaluates gradient of f.             

    Input Parameters:
.   tao  - the TAO_SOLVER context
.   X    - input vector
.   ptr  - optional user-defined context, as set by TaoAppSetConstraintRoutine()
    
    Output Parameters:
.   F - vector containing the newly evaluated gradient
*/
int FormConstraints(TAO_APPLICATION tao, Vec X, Vec F, void *ptr) 
{
  AppCtx *user = (AppCtx *) ptr;
  PetscScalar *x, *f;
  double *Vt1 = user->Vt1, *c = user->c, *d = user->d;
  double rate = user->rate;
  double dt = user->dt, ds = user->ds;
  PetscInt ms = user->ms;
  int    info;
  PetscInt i, xs,xm,gxs,gxm;
  Vec    localX,localF;
  double zero=0.0;

  info = DAGetLocalVector(user->da,&localX);CHKERRQ(info);
  info = DAGetLocalVector(user->da,&localF);CHKERRQ(info);
  info = DAGlobalToLocalBegin(user->da,X,INSERT_VALUES,localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(user->da,X,INSERT_VALUES,localX); CHKERRQ(info);
  info = DAGetCorners(user->da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(user->da,&gxs,PETSC_NULL,PETSC_NULL,&gxm,PETSC_NULL,PETSC_NULL); CHKERRQ(info);
  info = VecSet(F, zero);CHKERRQ(info);
  /* 
     The problem size is smaller than the discretization because of the
     two fixed elements (V(0,T) = E and V(Send,T) = 0.
  */

  /* Get pointers to the vector data */
  info = VecGetArray(localX, &x); CHKERRQ(info);
  info = VecGetArray(localF, &f); CHKERRQ(info);
  
  /* Left Boundary */
  if (gxs==0){ 
    f[0] = x[0]-user->strike;
  } else {
    f[0] = 0;
  }

  /* All points in the interior */
  /*  for (i=gxs+1;i<gxm-1;i++){ */
  for (i=1;i<gxm-1;i++){
    f[i] = (1.0/dt + rate)*x[i] - Vt1[i]/dt +
      (c[i]/(4*ds))*(x[i+1] - x[i-1] + Vt1[i+1] - Vt1[i-1]) +
      (d[i]/(2*ds*ds))*(x[i+1] -2*x[i] + x[i-1] + 
			Vt1[i+1] - 2*Vt1[i] + Vt1[i-1]);
  }

  /* Right boundary */
  if (gxs+gxm==ms){ 
    f[xm-1]=x[gxm-1];
  } else {
    f[xm-1]=0;
  }

  /* Restore vectors */
  info = VecRestoreArray(localX, &x); CHKERRQ(info);
  info = VecRestoreArray(localF, &f); CHKERRQ(info);
  info = DALocalToGlobalBegin(user->da,localF,F); CHKERRQ(info);
  info = DALocalToGlobalEnd(user->da,localF,F); CHKERRQ(info);
  info = DARestoreLocalVector(user->da,&localX); CHKERRQ(info);
  info = DARestoreLocalVector(user->da,&localF); CHKERRQ(info);
  info = PetscLogFlops(24*(gxm-2)); CHKERRQ(info);
  /*
  info=VecView(F,PETSC_VIEWER_STDOUT_WORLD);
  */
  return 0;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
/*
   FormJacobian - Evaluates Jacobian matrix.

   Input Parameters:
.  tao  - the TAO_SOLVER context
.  X    - input vector
.  ptr  - optional user-defined context, as set by TaoSetJacobian()

   Output Parameters:
.  J    - Jacobian matrix
*/
int FormJacobian(TAO_APPLICATION tao, Vec X, Mat *tJ, Mat *tJPre, MatStructure *flag, void *ptr)
{ 
  AppCtx *user = (AppCtx *) ptr;
  Mat J = *tJ;
  double *c = user->c, *d = user->d;
  double rate = user->rate;
  double dt = user->dt, ds = user->ds;
  PetscInt ms = user->ms;
  PetscScalar val[3];
  int info;
  PetscInt col[3];
  PetscInt i;
  PetscInt gxs,gxm;
  PetscTruth assembled;

  /* Set various matrix options */
  *flag=SAME_NONZERO_PATTERN;
  info = MatSetOption(J,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE); CHKERRQ(info);
  info = MatAssembled(J,&assembled); CHKERRQ(info);
  if (assembled){info = MatZeroEntries(J);  CHKERRQ(info);}


  info = DAGetGhostCorners(user->da,&gxs,PETSC_NULL,PETSC_NULL,&gxm,PETSC_NULL,PETSC_NULL); CHKERRQ(info);

  if (gxs==0){
    i = 0;
    col[0] = 0;
    val[0]=1.0;
    info = MatSetValues(J,1,&i,1,col,val,INSERT_VALUES); CHKERRQ(info);
  }
  for (i=1; i < gxm-1; i++) {
    col[0] = gxs + i - 1;
    col[1] = gxs + i;
    col[2] = gxs + i + 1;
    val[0] = -c[i]/(4*ds) + d[i]/(2*ds*ds);
    val[1] = 1.0/dt + rate - d[i]/(ds*ds);
    val[2] =  c[i]/(4*ds) + d[i]/(2*ds*ds);
    info = MatSetValues(J,1,&col[1],3,col,val,INSERT_VALUES); CHKERRQ(info);
  }
  if (gxs+gxm==ms){
    i = ms-1;
    col[0] = i;
    val[0]=1.0;
    info = MatSetValues(J,1,&i,1,col,val,INSERT_VALUES); CHKERRQ(info);
  }

  /* Assemble the Jacobian matrix */
  info = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = PetscLogFlops(18*(gxm)+5); CHKERRQ(info);
  return 0;
}

