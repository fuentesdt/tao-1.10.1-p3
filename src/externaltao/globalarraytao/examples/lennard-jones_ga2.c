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


/* Program usage: mpirun -machinefile m -np 1 lennard-jones_ga2 [-help] [all TAO options] */
/* Program usage: mpirun -machinefile m -np 4 lennard-jones_ga2 [-help] [all TAO options] */

/* Include MPI header files */
#include <mpi.h>

/*  Include tao header files so we can use TAO  */
#include "tao.h"
#include "taoapp_ga.h"

/* Include ga header files so that we can use GA */
#include "ga.h"
#include "taovec_ga.h"
#include "macdecls.h"


/*some standard header files */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>		//for malloc function
#include <math.h>


#define NDIM 2
#define NATOMS 1000
#define BLOCKSIZE 500
#define MAX_BLOCKS 25600

static char help[] = "This example demonstrates use of the TAO package to \n\
solve an unconstrained minimization problem on multiple processors.  We \n\
minimize the extended Lennard-Jones function.\n";


/*T 
   Concepts: TAO - Solving an unconstrained minimization problem
   Routines: TaoInitialize(); TaoFinalize(); TaoSetFromOptions();
   Routines: TaoGAApplicationCreate(); TaoSetApplication();
   Routines: TaoCreate(); TaoSetGAFunctionGradient(); 
   Routines: TaoSetGAHessian(); TaoSetGAInitialVector(); 
   Routines: TaoSolve(); TaoDestroy(); TaoApplicationDestroy();
   Routines: TaoGetTerminationReason();
   Processors: >=1 
T*/


/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines that evaluate the function,
   gradient, and hessian.
*/
typedef struct {
  int x; 
  int y;
} topo_t;
typedef struct
{
  int ndim;			/* dimension n = NDIM*NATOMS */
  int natoms;                   /* number of atoms */
  int me;                       /* proc id  */
  int nproc;                    /* number of processes */
  int BlockSize;
  int nBlocks;
  int memHandle;               /* memory handle for MA allocated work space */
  int TopoHandle;              /* memory handle for allocated btopo */
  double *x1,*x2,*grad;         /* work arrays */
  topo_t btopo[MAX_BLOCKS];                /* block topography array  */
  int atomicTask;               /* Atomic Task: 1-d array whose size=1 */
} AppCtx;


/* -------------- User-defined routines ---------- */

int FormFunctionGradient (TAO_GA_APPLICATION gaapp, GAVec ga_X, double *f, 
			  GAVec ga_G, void *ptr);
int InitializeVariables(GAVec ga_x, AppCtx *user);
int SetupBlocks(AppCtx *user);
int getBlock(GAVec ga_x, int taskId, AppCtx *user);

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
  int heap = 400000, stack = 400000;
  MPI_Init (&argc, &argv);	/* initialize MPI */
  GA_Initialize ();		/* initialize GA */
  user.me = GA_Nodeid ();
  user.nproc = GA_Nnodes ();
  startTime = MPI_Wtime();
  
  if (user.me == 0) {
    if (GA_Uses_fapi ())
      GA_Error ("Program runs with C array API only", 0);
    printf ("Using %ld processes\n", (long) user.nproc);
    fflush (stdout);
  }
  heap /= user.nproc;
  stack /= user.nproc;
  if (!MA_init (MT_F_DBL, stack, heap))
    GA_Error ("MA_init failed", stack + heap);	/* initialize memory allocator */
  
  /* Initialize TAO */
  TaoInitialize (&argc, &argv, (char *) 0, help);

  /* Initialize problem parameters */
  user.ndim = NDIM;
  user.natoms = NATOMS;
  user.BlockSize = BLOCKSIZE;


  /* Allocate vectors for the solution and gradient */
  int dims[2];
  dims[0] = user.ndim*user.natoms;
  ga_x = NGA_Create (C_DBL, 1, dims, "GA_X", NULL);
  if (!ga_x) GA_Error ("lennard-jones.main::NGA_Create ga_x", ga_x);

  /* Set up structures for data distribution */
  info = SetupBlocks(&user); CHKERRQ(info);


  /* The TAO code begins here */
  /* Create TAO solver with desired solution method */
  info = TaoCreate (MPI_COMM_WORLD, "tao_lmvm", &tao); CHKERRQ(info);
  info = TaoGAApplicationCreate (MPI_COMM_WORLD, &taoapp); CHKERRQ(info);

  /* Set the initial solution */
  info = InitializeVariables(ga_x, &user); CHKERRQ(info);
  info = TaoGAAppSetInitialSolutionVec(taoapp, ga_x); CHKERRQ(info);

  /* Set routines for function, gradient */
  info = TaoGAAppSetObjectiveAndGradientRoutine (taoapp, FormFunctionGradient, 
					      (void *) &user); CHKERRQ(info);

  /* Check for TAO command line options */
  info = TaoSetFromOptions (tao); CHKERRQ(info);


  /* SOLVE THE APPLICATION */
  info = TaoSolveGAApplication (taoapp, tao); CHKERRQ(info);

  /*  To View TAO solver information use */
  info = TaoView(tao); CHKERRQ(info);

  /* Get termination information */
  info = TaoGetTerminationReason (tao, &reason);
  if(info) GA_Error("lennard-jones.main.TaoGetTerminationReason",info);
  if (user.me == 0) {
    if (reason <= 0)
      printf("Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
    
    printf("WALL TIME TAKEN = %lf\n", MPI_Wtime()-startTime);
    /*output the solutions */ 

    printf ("The solution is :\n");
  }
  GA_Print (ga_x);




  /* Free TAO data structures */
  info = TaoDestroy (tao); CHKERRQ(info);
  info = TaoGAAppDestroy (taoapp); CHKERRQ(info);

  /* Free GA data structures */
  GA_Destroy (ga_x);
  if (!MA_pop_stack(user.memHandle)) 
    ga_error("Main::MA_pop_stack for memHandle failed",0);

  /* Finalize TAO, GA, and MPI */
  TaoFinalize ();
  GA_Terminate ();
  MPI_Finalize ();

  return 0;
}



#define FUNCTION_GRAD_2D(is, js) \
      xx = x_j[2*j] - x_i[2*i];           \
      yy = x_j[2*j+1] - x_i[2*i+1];       \
      rij = xx*xx + yy*yy;                \
      temp = 1.0/rij/rij/rij;             \
      *f += temp*(temp-2.0);              \
      temp *= 12.0*(temp-1.0)/rij;        \
      g[js + 2*j] -= xx*temp;             \
      g[js + 2*j + 1] -= yy*temp;         \
      g[is+2*i] += xx*temp;               \
      g[is + 2*i +1] += yy*temp;          

int LJFG_2D(int taskId,double *f,AppCtx *user)
{

  int b_x, b_y; /* block topology */
  int i, j, start_i, start_j;
  double xx, yy, rij, temp;
  double *g, *x_i, *x_j;

  b_x = user->btopo[taskId].x;
  b_y = user->btopo[taskId].y;



  start_i = user->BlockSize * b_x * 2;
  start_j = user->BlockSize * b_y * 2;
  g = user->grad;
  x_i = user->x1;
  x_j = user->x2;

  if (b_x == b_y)  /* compute lower triangular matrix only */
    for (i=1; i<user->BlockSize; i++) 
      for (j=0; j<i; j++) {
	FUNCTION_GRAD_2D(start_i, start_j);
      }
  else if (b_x > b_y)  /* compute right half of the block */
    for (i=0; i<user->BlockSize; i++)
      for (j=user->BlockSize/2; j<user->BlockSize; j++) {
	FUNCTION_GRAD_2D(start_i, start_j);
      }
  else /* Compute upper half of the block */
    for (i=0; i<user->BlockSize/2; i++) 
      for (j=0; j<user->BlockSize; j++) {
	FUNCTION_GRAD_2D(start_i, start_j);
      }
  return 0;

}



#define FUNCTION_GRAD_3D(is, js) \
      xx = x_j[3*j] - x_i[3*i];           \
      yy = x_j[3*j+1] - x_i[3*i+1];       \
      zz = x_j[3*j+2] - x_i[3*i+2];       \
      rij = xx*xx + yy*yy + zz*zz;        \
      temp = 1.0/rij/rij/rij;             \
      *f += temp*(temp-2.0);              \
      temp *= 12.0*(temp-1.0)/rij;        \
      g[js + 3*j] -= xx*temp;             \
      g[js + 3*j + 1] -= yy*temp;         \
      g[js + 3*j + 2] -= zz*temp;         \
      g[is + 3*i] += xx*temp;               \
      g[is + 3*i +1] += yy*temp;          \
      g[is + 3*i +2] += zz*temp;          

int LJFG_3D(int taskId,double *f,AppCtx *user)
{

  int b_x, b_y; /* block topology */
  int i, j, start_i, start_j;
  double xx, yy, zz, rij, temp;
  double *g, *x_i, *x_j;

  b_x = user->btopo[taskId].x;
  b_y = user->btopo[taskId].y;
  
  start_i = user->BlockSize * b_x * 3;
  start_j = user->BlockSize * b_y * 3;

  g = user->grad;
  x_i = user->x1;
  x_j = user->x2;

  if (b_x == b_y)  /* compute lower triangular matrix only */
    for (i=1; i<user->BlockSize; i++) 
      for (j=0; j<i; j++) {
	FUNCTION_GRAD_3D(start_i, start_j);
      }
  else if (b_x > b_y)  /* compute right half of the block */
    for (i=0; i<user->BlockSize; i++)
      for (j=user->BlockSize/2; j<user->BlockSize; j++) {
	FUNCTION_GRAD_3D(start_i, start_j);
      }
  else /* Compute upper half of the block */
    for (i=0; i<user->BlockSize/2; i++) 
      for (j=0; j<user->BlockSize; j++) {
	FUNCTION_GRAD_3D(start_i, start_j);
      }
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
  AppCtx *user = (AppCtx *) ptr;
  int lo, hi;		
  int taskId=user->me;          //Which task am I running
  int i;
  int zero = 0;

  /* reset atomicTask to nproc */
  if (user->me == 0) 
    NGA_Put(user->atomicTask, &zero, &zero, &user->nproc, &zero);


  for (i=0;i<user->natoms*user->ndim; i++)
    user->grad[i] = 0.0;
  *f = 0.0;

  while (taskId < user->nBlocks) {
    getBlock(ga_X, taskId, user);
    if (user->ndim == 2)
      LJFG_2D(taskId,f, user);
    else
      LJFG_3D(taskId,f, user);

    /* Get next block */
    taskId += user->nproc;
    //NGA_Read_inc(user->atomicTask, &zero, 1); 
  }


  /* gather function */
  GA_Dgop(f, 1, "+");

  /* gather gradient */
  GA_Dgop(user->grad, user->natoms*user->ndim, "+");
  NGA_Distribution(ga_G, user->me, &lo, &hi);
  NGA_Put(ga_G, &lo, &hi, user->grad+lo, &hi);

  GA_Sync();
  
  return 0;
}


int InitializeVariables(GAVec ga_x, AppCtx *user) {
  double *x;
  int lo, hi, n, handle;
  double xx,yy,zz;
  int i,j,k, il, jl, ctr, icrtn, ileft, isqrtn;

  // This method is not parallelized
  if (user->me != 0) {
    GA_Sync();
    return 0;
  }

  n = user->ndim * user->natoms + 1;
  if (MA_push_stack(C_DBL, n, "InitializeVariables buffer", &handle))
    MA_get_pointer(handle, &x);
  else 
    ga_error("ma_alloc_get failed", n);

  lo = 0;
  hi = n-2;

  if (user->ndim == 2) {
    isqrtn = int(sqrt(user->natoms));
    ileft = user->natoms - isqrtn * isqrtn;
    xx = yy = 0.0;
    for (j=0; j <= isqrtn + ileft/isqrtn; j++) {
      for (i=0; i<TaoMin(isqrtn, user->natoms - j * isqrtn); i++) {
	ctr = j * isqrtn + i;
	x[2*ctr] = xx;
	x[2*ctr + 1] = yy;
	xx += 1.0;
      }
      yy += 1.0;
      xx = 0.0;
    }
  } else if (user->ndim == 3) {
    icrtn = (int)pow((int)(user->natoms + 0.5),1.0/3.0);
    ileft = user->natoms - icrtn*icrtn*icrtn;
    xx = yy = zz = 0.0;
    for (k=0; k<=icrtn + ileft/(icrtn*icrtn); k++) {
      jl = TaoMin(icrtn, (user->natoms - k*icrtn*icrtn)/icrtn + 1);
      for (j=0; j<jl; j++) {
	il = TaoMin(icrtn, user->natoms - k*icrtn*icrtn - j*icrtn);
	for (i = 0; i<il; i++) {
	  ctr = k*icrtn*icrtn + j*icrtn + i;
	  x[3*ctr] = xx;
	  x[3*ctr + 1] = yy;
	  x[3*ctr + 2] = zz;
	  xx += 1.0;
	}
	yy += 1.0;
	xx = 0.0;
      }
      zz += 1.0;
      yy = 0.0;
    }
  }

  // Distribute the array
  NGA_Put(ga_x, &lo, &hi, x, &hi);

  if (!MA_pop_stack(handle)) 
    ga_error("InitializeVariables:MA_pop_stack failed",0);

  GA_Sync();

  return 0;
}


/**
 * Block Topology (of Force Matrix): 
 * Say for example: If there are 4 block and 100 atoms, the size of 
 * the force matrix is 100x100 and each block size is 50x50. 
 * 
 *  -----------
 * |     |     |
 * | 0,0 | 0,1 |
 *  -----------
 * | 1,0 | 1,1 |
 * |     |     |
 *  -----------
 */
int SetupBlocks(AppCtx *user)
{
  int i,j,k=0;
  int n;
  int zero = 0;
  int x_space, g_space;

  if (user->natoms % user->BlockSize) {
    GA_Error("Number of atoms should be a multiple of block size. Choose a different block size.", 0L);
  }

  n = user->natoms / user->BlockSize;
  user->nBlocks = n*n;

  if (user->nBlocks > MAX_BLOCKS) 
    GA_Error("Number of blocks is greater that MAX_BLOCKS: Solution is either to increase the defined MAX_BLOCKS or increase your block size",0L);

  if (user->nBlocks < user->nproc)
    GA_Error("Number of blocks should be greater than or equal to the number of processors",0L);

  
  for (i=0;i<n;i++)
    for (j=0;j<n;j++,k++) {
      user->btopo[k].x = i;
      user->btopo[k].y = j;
    }
  
  /* Create task array */
  n = 1;
  user->atomicTask = NGA_Create(C_INT, 1, &n, "Atomic Task", NULL);
  if (!user->atomicTask)
    GA_Error("NGA_Create failed for Atomic Task",0);

  if (user->me == 0) 
    NGA_Put(user->atomicTask, &zero, &zero, &user->nproc, &zero);
  
  
  /* space for x values from two processors */
  x_space = 2 * user->BlockSize * user->ndim;
  /* space for ALL gradient value */
  g_space = user->natoms * user->ndim; 
             

  if (MA_push_stack(C_DBL, x_space + g_space+3, "GA LJ bufs", &user->memHandle))
    MA_get_pointer(user->memHandle, &user->x1);
  else
    GA_Error("ma_alloc_get failed",x_space + g_space);
  
  user->x2  = user->x1 + x_space/2 + 1;
  user->grad = user->x2 + x_space/2 + 1;
  GA_Sync();
  return 0;
}



int getBlock(GAVec g_x, int taskId, AppCtx *user) 
{
  int lo, hi;
  int size;

  size = user->BlockSize * user->ndim;

  /* get the coordinates of the atoms in the corresponding rows in the block */
  lo = user->btopo[taskId].x*size;
  hi = lo + size -1;
  NGA_Get(g_x, &lo, &hi, user->x1, &hi);

  /* get the coordinates of the atoms in the corresponding cols in the block */
  lo = user->btopo[taskId].y*size;
  hi = lo + size - 1;
  NGA_Get(g_x, &lo, &hi, user->x2, &hi);
  return 0;
}
