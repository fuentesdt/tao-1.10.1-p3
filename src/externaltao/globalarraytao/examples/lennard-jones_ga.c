/*File: lennard-jones_ga.c
  Purpose: to test Lennard-Jones problem using Tao/GA vectors using multiple processors.
  Date created: August 16, 2002
*/
     

/**************************************************************

Author: Limin Zhang, Ph.D.
        Mathematics Department
        Columbia Basin College
        Pasco, WA 99301
        Limin.Zhang@cbc2.org

Mentor: Jarek Naplocha, Ph.D.
        Environmental Molecular Science Laboratory
        Pacific Northwest National Laboratory
        Richland, WA 99352

Date: 4/22/2002

Purpose:
      to design and implement an application for  between TAO and
      global arrays.
**************************************************************/

/* Include MPI header files */
#include <mpi.h>

/* Program usage: mpirun -machinefile m -np 1 lennard-jones_ga [-help] [all TAO options] */
/* Program usage: mpirun -machinefile m -np 4 lennard-jones_ga [-help] [all TAO options] */

/*  Include tao header files so we can use TAO  */
#include "tao.h"
#include "taoapp_ga.h"

/* Include ga header files so that we can use GA */
#include "ga.h"
#include "taovec_ga.h"
#include "macdecls.h"


/*some standard header files */
#include <stdio.h>
#include <math.h>

#define NDIM 2
#define NATOMS 16

static char help[] = "This example demonstrates use of the TAO package to \n\
solve an unconstrained minimization problem on multiple processors.  We \n\
minimize the extended Lennard-Jones function.";

/*T 
   Concepts: TAO - Solving an unconstrained minimization problem
   Routines: TaoInitialize(); TaoFinalize(); TaoSetFromOptions();
   Routines: TaoGAApplicationCreate(); TaoGASetInitialSolutionVec();
   Routines: TaoCreate(); TaoGASetObjectiveAndGradientRoutine();
   Routines: TaoSetGAHessian(); TaoSetGAInitialVector(); 
   Routines: TaoSolveGAApplication(); TaoDestroy(); TaoGAAppDestroy();
   Routines: TaoGetTerminationReason(); TaoView();
   Processors: 1 
T*/


/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines that evaluate the function,
   gradient, and hessian.
*/
typedef struct
{
  int n;			/* dimension n = NDIM*NATOMS */
  int ndim;
  int natoms;
  int memHandle;
}
AppCtx;


/* -------------- User-defined routines ---------- */

int FormFunctionGradient (TAO_GA_APPLICATION gaapp, GAVec ga_X, double *f, GAVec ga_G, void *ptr);
int InitializeVariables(GAVec ga_X, AppCtx *user); 


int main (int argc, char **argv)
{
  double startTime;
  int info;			/* used to check for functions returning nonzeros */
  GAVec ga_x;        		/* solution vector */
  TAO_SOLVER tao;		/* TAO_SOLVER solver context */
  TAO_GA_APPLICATION taoapp;	/* TAO application context */
  TaoTerminateReason reason;
  AppCtx user;			/* user-defined application context */

  /*initialize GA and MPI */
  int heap = 4000, stack = 4000;
  MPI_Init (&argc, &argv);	/* initialize MPI */
  GA_Initialize ();		/* initialize GA */
  if (!MA_init(MT_F_DBL, stack, heap))
    GA_Error((char*)"MA_init failed", stack+heap);

  /* Initialize TAO */
  TaoInitialize (&argc, &argv, (char *) 0, help);

  startTime = MPI_Wtime();

  /* Initialize problem parameters */
  user.natoms = NATOMS;
  user.ndim = NDIM;
  user.n = user.natoms*user.ndim;

  /* Create working space */
  if (MA_push_stack(C_DBL, 2*user.n, "Vector buffers", &user.memHandle) == MA_FALSE)
    GA_Error((char*)"MAIN::ma_alloc_get failed",2*user.n);

  /* Allocate Global Array vector for the solution */
  int dims[2];
  dims[0] = user.n;
  ga_x = NGA_Create (C_DBL, 1, dims, (char*)"GA_X", NULL);
  if (!ga_x) GA_Error ((char*)"lennard-jones.main::NGA_Create ga_x", ga_x);

  /* The TAO code begins here */
  /* Create TAO solver with desired solution method */
  info = TaoCreate (MPI_COMM_WORLD, "tao_cg", &tao); CHKERRQ(info);
  info = TaoGAApplicationCreate (MPI_COMM_WORLD, &taoapp); CHKERRQ(info);

  /* Set initial vector */
  info = InitializeVariables(ga_x, &user); CHKERRQ(info);
  info = TaoGAAppSetInitialSolutionVec(taoapp, ga_x); CHKERRQ(info);

  /* Set routines for function, gradient */
  info = TaoGAAppSetObjectiveAndGradientRoutine (taoapp, FormFunctionGradient, (void *) &user); 
  CHKERRQ(info);

  /* Check for TAO command line options */
  info = TaoSetFromOptions (tao); CHKERRQ(info);

  /* SOLVE THE APPLICATION */
  info = TaoSolveGAApplication (taoapp, tao); CHKERRQ(info);

  /*  To View TAO solver information use */
  info = TaoView(tao); CHKERRQ(info);

  /* Get termination information */
  info = TaoGetTerminationReason (tao, &reason); CHKERRQ(info);

  if (reason <= 0)
    printf("Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");

  printf("TIME TAKEN = %lf\n", MPI_Wtime()-startTime);

  /*output the solutions */ 
  printf ("The solution is :\n");
  GA_Print (ga_x);

  /* Free TAO data structures */
  info = TaoDestroy (tao); CHKERRQ(info);
  info = TaoGAAppDestroy (taoapp); CHKERRQ(info);

  /* Free GA data structures */
  GA_Destroy (ga_x);
  if (!MA_pop_stack(user.memHandle))
    GA_Error((char*)"Main::MA_pop_stack failed",0);

  /* Finalize TAO, GA, and MPI */
  TaoFinalize ();
  GA_Terminate ();
  MPI_Finalize ();

  return 0;
}

/* -------------------------------------------------------------------- */
/*  
    FormFunctionGradient - Evaluates the function, f(X), and gradient, G(X). 

    Input Parameters:
.   tao  - the TAO_SOLVER context
.   ga_X    - input vector
.   ptr  - optional user-defined context, as set by TaoSetFunctionGradient()
    
    Output Parameters:
.   ga_G - vector containing the newly evaluated gradient
.   f - function value

    Note:
    Some optimization methods ask for the function and the gradient evaluation
    at the same time.  Evaluating both at once may be more efficient that
    evaluating each separately. 
*/


int FormFunctionGradient (TAO_GA_APPLICATION gaapp, GAVec ga_X, double *f, GAVec ga_G, void *ptr)
{
  int lo, hi;			//the global coordinates
  AppCtx *user = (AppCtx *) ptr;
  int i,j;
  double *g, *x;
  double xx,yy,zz,temp,rij;
  

  MA_get_pointer(user->memHandle, &x);
  g = x + user->ndim*user->natoms;
    
  lo=0;
  hi=user->n-1; /* range of array indices */

  NGA_Get(ga_X, &lo, &hi, x, &hi);

  *f = 0;
  for (i=0; i < user->n; i++)
    g[i] = 0.0;

  if (user->ndim == 2) {
    for (j=1; j < user->natoms; j++) {
      for (i=0; i<j; i++) {
	xx = x[2*j] - x[2*i];
	yy = x[2*j+1] - x[2*i+1];
	rij = xx*xx + yy*yy;
	temp = 1.0/rij/rij/rij;
	*f += temp*(temp-2.0);
	temp *= 12.0*(temp-1.0)/rij;
	g[2*j] -= xx*temp;
	g[2*j+1] -= yy*temp;
	g[2*i] += xx*temp;
	g[2*i+1] += yy*temp;
      }
    }
  } else if (user->ndim == 3) {
    for (j=1; j < user->natoms; j++) {
      for (i=0; i < j; i++) {
	xx = x[3*j] - x[3*i];
	yy = x[3*j+1] - x[3*i+1];
	zz = x[3*j+2] - x[3*i+2];
	rij = xx*xx + yy*yy + zz*zz;
	temp = 1.0/rij/rij/rij;
	*f += temp*(temp-2.0);
	temp *= 12.0*(temp-1.0)/rij;
	g[3*j] -= xx*temp;
	g[3*j+1] -= yy*temp;
	g[3*j+2] -= zz*temp;
	g[3*i] += xx*temp;
	g[3*i+1] += yy*temp;
	g[3*i+2] += zz*temp;
      }
    }
  }
      
  NGA_Put(ga_G, &lo, &hi, g, &hi);
  return 0;
}

int InitializeVariables(GAVec ga_X, AppCtx *user) 
{
  double *x;
  double xx, yy, zz;
  int isqrtn, icrtn, left, i, j, k, il, jl, ctr;
  int lo, hi;


  lo = 0; hi = user->n - 1; /* range of array indices */
  MA_get_pointer(user->memHandle, &x);

  NGA_Get (ga_X, &lo, &hi, x, &hi);

  if (user->ndim == 2) {
    isqrtn = (int) sqrt( (double) user->natoms);
    left = user->natoms - isqrtn * isqrtn;
    xx = 0.0;
    yy = 0.0;
    for (j=0; j<=isqrtn + left/isqrtn; j++) {
      for (i=0; i < TaoMin(isqrtn, user->natoms - j*isqrtn); i++) {
	ctr = j*isqrtn + i;
	x[2*ctr] = xx;
	x[2*ctr+1] = yy;
	xx += 1.0;
      }
      yy += 1.0;
      xx = 0.0;
    }	     
  }
  else if (user->ndim == 3) {
    icrtn = (int) pow((user->natoms + 0.5),1.0/3.0);
    left = user->natoms - icrtn * icrtn * icrtn;
    xx = yy = zz = 0.0;
    for (k=0; k <= icrtn + left; k++) {
      jl = TaoMin(icrtn, (user->natoms - k*icrtn*icrtn)/icrtn+1);
      for (j=0; j<jl; j++) {
	il = TaoMin(icrtn, user->natoms - k*icrtn*icrtn - j*icrtn);
	for (i=0; i<il; i++) {
	  ctr = k*icrtn*icrtn + j*icrtn + i;
	  x[3*ctr] = xx;
	  x[3*ctr+1] = yy;
	  x[3*ctr+2] = zz;
	  xx += 1.0;
	}
	yy += 1.0;
	xx = 0.0;
      }
      zz += 1.0;
      yy = 0.0;
    }
  }

  NGA_Put(ga_X, &lo, &hi, x, &hi);
  return 0;
}

