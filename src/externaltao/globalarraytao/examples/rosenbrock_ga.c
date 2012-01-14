/*File: rosenbrock1.c

  Purpose: 
		   it is the GA version of rosenbrock1.c to test the class TaoGAApplication. 
  Acknowledgement: 
		   the code was originally developed by Dr.Steve Bensor and others at ANL and modified by Dr. Limin Zhang
  		   and Dr. Jarek Nieplochai from PNNL.
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


/*$Id: s.rosenbrock1.c 1.41 02/02/14 14:22:13-06:00 sarich@holmes.mcs.anl.gov $*/

/* Program usage: mpirun -machinefile m -np 1 rosenbrock1 [-help] [all TAO options] */
/* Program usage: mpirun -machinefile m -np 4 rosenbrock1 [-help] [all TAO options] */
/*MPI header*/
#include <mpi.h>

/*  Include "tao.h" so we can use TAO solvers.  */
#include "tao.h"		//for tao initialization

/*include "ga.h" so that we can use GA*/
#include "ga.h"
#include "taovec_ga.h"
#include "taoapp_ga.h"
#include "macdecls.h"



/*some standard C header files */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>		//for malloc function
#include <math.h>


static char help[] = "This example demonstrates use of the TAO package to \n\
solve an unconstrained minimization problem on a single processor.  We \n\
minimize the extended Rosenbrock function: \n\
   sum_{i=0}^{n/2-1} ( alpha*(x_{2i+1}-x_{2i}^2)^2 + (1-x_{2i})^2 ) \n";

/*T 
   Concepts: TAO - Solving an unconstrained minimization problem
   Routines: TaoInitialize(); TaoFinalize(); TaoSetFromOptions();
   Routines: TaoGAApplicationCreate(); TaoSetApplication();
   Routines: TaoCreate(); TaoSetGAFunctionGradient(); 
   Routines: TaoSetGAHessian(); TaoSetGAInitialVector(); 
   Routines: TaoSolve(); TaoDestroy(); TaoApplicationDestroy();
   Routines: TaoGetTerminationReason();
   Processors: 1
T*/


/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines that evaluate the function,
   gradient, and hessian.
*/
typedef struct
{
  int n;			/* dimension */
  double alpha;			/* condition parameter */
}
AppCtx;

/* -------------- User-defined routines ---------- */
int FormFunctionGradient (TAO_GA_APPLICATION, GAVec, double *, GAVec, void *);
int SetVariableBounds (TAO_GA_APPLICATION, GAVec, GAVec, void*);

int main (int argc, char **argv)
{
  int info;			/* used to check for functions returning nonzeros */
  GAVec ga_x;			/* solution, gradient vectors */
  TAO_SOLVER tao;		/* TAO_SOLVER solver context */
  TAO_GA_APPLICATION taoapp;	/* TAO application context */
  TaoTerminateReason reason;
  AppCtx user;			/* user-defined application context */

  /*initialize GA and MPI */
  int heap = 20000, stack = 20000;
  int me, nproc;
  MPI_Init (&argc, &argv);	/* initialize MPI */
  GA_Initialize ();		/* initialize GA */
  me = GA_Nodeid ();
  nproc = GA_Nnodes ();
  if (me == 0) {
    if (GA_Uses_fapi ())
      GA_Error((char*)"Program runs with C array API only", 0);
    printf ("Using %ld processes\n", (long) nproc);
    fflush (stdout);
  }
  heap /= nproc;
  stack /= nproc;
  if (!MA_init (MT_F_DBL, stack, heap))
    GA_Error((char*)"MA_init failed", stack + heap);	/* initialize memory allocator */

  /* Initialize TAO */
  TaoInitialize (&argc, &argv, (char *) 0, help);

  /* Initialize problem parameters */
  user.n = 2;
  user.alpha = 99.0;

  /* Allocate vectors for the solution and gradient */
  int ndim = 1, dims[2], type = C_DBL;
  dims[0] = user.n;
  ga_x = NGA_Create (type, ndim, dims, (char*)"GA_X", NULL);
  if (!ga_x) GA_Error((char*)"rosenbrock1.main::NGA_Create ga_x", ga_x);
  /*  ga_g = GA_Duplicate (ga_x, "GA_G");
      if (!ga_g) GA_Error((char*)"rosenbrock1.main::NGA_Create ga_g", ga_g); */

  /* The TAO code begins here */
  /* Create TAO solver with desired solution method */
  info = TaoCreate (MPI_COMM_WORLD, "tao_blmvm", &tao); CHKERRQ(info);
  info = TaoGAApplicationCreate (MPI_COMM_WORLD, &taoapp); CHKERRQ(info);

  info = TaoGAAppSetInitialSolutionVec(taoapp, ga_x); CHKERRQ(info);

  info = TaoGAAppSetVariableBoundsRoutine(taoapp, SetVariableBounds, (void*)&user); CHKERRQ(info);

  /* Set routines for function, gradient */
  info = TaoGAAppSetObjectiveAndGradientRoutine (taoapp, FormFunctionGradient, (void *) &user);	 CHKERRQ(info);



  /* Check for TAO command line options */
  info = TaoSetFromOptions (tao); CHKERRQ(info);

  /* SOLVE THE APPLICATION */
  info = TaoSolveGAApplication (taoapp,tao); CHKERRQ(info);


  /*  To View TAO solver information use */
  info = TaoView(tao); CHKERRQ(info);

  

  /* Get termination information */
  info = TaoGetTerminationReason (tao, &reason); CHKERRQ(info);

  if (reason <= 0)
    printf
      ("Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");

  /*output the solutions */ 
  if (me == 0)
    printf ("The solution is :\n");
  GA_Print (ga_x);


  /* Free TAO data structures */
  info = TaoDestroy (tao); CHKERRQ(info);
  info = TaoGAAppDestroy (taoapp); CHKERRQ(info);


  /* Free GA data structures */
  GA_Destroy (ga_x);

  /* Finalize TAO, GA, and MPI */
  TaoFinalize ();
  GA_Terminate ();
  MPI_Finalize ();

  if (me == 0)
    printf ("Normal exit...\n");

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


int
FormFunctionGradient (TAO_GA_APPLICATION gaapp, GAVec ga_X, double *f, GAVec ga_G, void *ptr)
{
  AppCtx *user = (AppCtx *) ptr;
  int i, nn = user->n / 2;
  double ff = 0., t1, t2, alpha = user->alpha;
  double glocal[2],xlocal[2];

  int me = GA_Nodeid ();

  GA_Sync();

  if(me ==0){
    
     /* given the size of the problem we will execute this part on process 0 only */

     int lo=0,hi=1; /* range of array indices */
     NGA_Get (ga_X, &lo, &hi, xlocal, &hi);
     /* Compute G(X) */

     for (i=0; i<nn; i++){
       t1 = xlocal[2*i+1]-xlocal[2*i]*xlocal[2*i]; t2= 1-xlocal[2*i];
       ff += alpha*t1*t1 + t2*t2;
       glocal[2*i] = -4*alpha*t1*xlocal[2*i]-2.0*t2;
       glocal[2*i+1] = 2*alpha*t1;
     }
   
     GA_Init_fence();
     NGA_Put(ga_G, &lo, &hi, glocal, &hi);
     GA_Fence();
  }   

  /* broadcast the result to everybody */
  GA_Brdcst(&ff, sizeof(double), 0);

  *f = ff;
  return 0;
}


int SetVariableBounds(TAO_GA_APPLICATION gaapp, GAVec ga_xl, GAVec ga_xu, void *ctx)
{
  double lb = 0.0, ub = 2.0;
  
  GA_Fill(ga_xl, (void*)&lb);
  GA_Fill(ga_xu, (void*)&ub);
  return 0;
}
