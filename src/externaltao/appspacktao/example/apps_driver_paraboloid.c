/* apps_driver.c
 *
 *
 *  Find the minimizer of f(x,y) = (1.5-x(1-y))^2 + (2.25 - x(1-y)^2)^2 + (2.625 - x(1-y^3))^3
 */

/* Or alternatively, find the minimizer of of f(x,y) = (x-4)^2 + (y-3)^2 + (z+7)^2 */


#include <math.h>
#include "tao_solver.h"
#include "tao.h"
#include "src/petsctao/matrix/taomat_petsc.h"
#include "src/petsctao/vector/taovec_petsc.h"
#include "taovec.h"
#include "taomat.h"
#include "apps_solver.h"
#include "petsc.h"



extern int TaoCreate_APPS(TAO_SOLVER);
static const char help[] = "No help available";

//extern "C" int TaoCreate_APPS(TAO_SOLVER);

typedef struct {
  int n;
  double alpha;
  double c;
  int argc;
  char **argv;
} AppCtx;

/* These are necessary user-defined functions */
int funcgrad(TAO_SOLVER solver, Vec X, double *f, Vec GG, void *ctx);
int hessian(TAO_SOLVER solver, Vec X, Mat *M, Mat *M, MatStructure *flg, void *ctx); 
//int funcgrad2(TAO_SOLVER solver, TaoVec *XX, double *f, TaoVec *GG, void *ctx);
//int hessian2(TAO_SOLVER solver, TaoVec *XX, TaoMat *MM, void *ctx); 
int my_monitor(TAO_SOLVER solver, void *mctx);

/* Defined in apps_solver.C */
extern int setargs(TAO_SOLVER solver, int argc, char *argv[]);

// rosenbrock1 functions
int FormFunctionGradient(TAO_SOLVER, TaoVec*, double*, TaoVec*, void*);
int FormHessian(TAO_SOLVER, TaoVec*, TaoMat*, void*);

//==========================================
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  int info;
  TaoTruth success; // flag
  PetscTruth flg;
  TAO_SOLVER tao;
  TAO_APPLICATION taoapp;
  Vec X,G;
  Mat H;
  Vec L,U;
  double f;
  double *x0;     // starting vector array
  int ind[2]={0,1};  // indices for starting vector
  PetscScalar sx, sy;   // starting vector coordinates
  Vec X0;    // starting Vec
  double *x,*l,*u;
  AppCtx user;
  

  TaoFunctionBegin;

  user.n = 2; user.alpha = 99;

  // Initialize
  info = PetscInitialize(&argc, &argv, (char *)0, help); CHKERRQ(info);
  info = TaoInitialize(&argc, &argv, (char *)0,help); CHKERRQ(info);


  // Add the apps method to the dynamic registry
  info = TaoRegisterAppspack(); CHKERRQ(info);
  //info = TaoRegisterDynamic("tao_apps","/home/sarich/testtao/tao/lib/libg_c++/linux/libappspacktao.so","TaoCreate_APPS",TaoCreate_APPS); CHKERRQ(info);

  sx = sy = 0;   // default
  info = PetscOptionsGetScalar((char*)0,"-sx",&sx,&flg); CHKERRQ(info);
  info = PetscOptionsGetScalar((char*)0,"-sy",&sy,&flg); CHKERRQ(info);
  
  
  
  // Create the vectors X, G, X0
  info = VecCreateSeq(PETSC_COMM_SELF,user.n, &X); CHKERRQ(info);

  info = VecDuplicate(X,&G); CHKERRQ(info);
  info = VecDuplicate(X,&X0); CHKERRQ(info);

  // Get a handle to the vector X0, set to [sx, sy], return the handle
  info = VecGetArray(X0, &x0); CHKERRQ(info);
  x0[0] = sx; x0[1] = sy; 
  info = VecRestoreArray(X0, &x0); CHKERRQ(info);

  // Create the hessian
  info = MatCreateSeqAIJ(PETSC_COMM_SELF,user.n, user.n, user.n, PETSC_NULL,&H); CHKERRQ(info);


  // Create the tao object
  info = TaoCreate(PETSC_COMM_WORLD,"tao_apps",&tao);CHKERRQ(info);

  info = TaoPetscApplicationCreate(PETSC_COMM_WORLD,&taoapp);

  // Set the functions 
  info = TaoSetPetscFunctionGradient(taoapp,X,G,funcgrad,(void*)&user); CHKERRQ(info);
  info = TaoSetPetscHessian(taoapp,H,H,hessian,(void*)&user); CHKERRQ(info);


  //info = TaoSetMonitor(tao,my_monitor,(void*)0,0); CHKERRQ(info);

  info = TaoSetApplication(tao,taoapp); CHKERRQ(info);
  // Use the PetscVec X0 as the Initial TaoVec XX0
  info = TaoSetPetscInitialVector(taoapp,X0); CHKERRQ(info);

  // Set bounds here
  info = VecDuplicate(X,&L); CHKERRQ(info);
  info = VecDuplicate(X,&U); CHKERRQ(info);

  info = VecGetArray(L,&l); CHKERRQ(info);
  l[0] = -4;
  l[1] = -1;
  info = VecRestoreArray(L,&l); CHKERRQ(info);
  
  info = VecGetArray(U,&u); CHKERRQ(info);
  u[0] = 1;
  u[1] = 1;
  info = VecRestoreArray(U,&u); CHKERRQ(info);


  //info = TaoSetVariableBounds(tao,LL,UU); CHKERRQ(info);
  /* Send argv and argc to tao */
  info = setargs(tao,argc, argv);

  // Actually solve the problem
  info = TaoSolve(tao); CHKERRQ(info);

  info = TaoView(tao); CHKERRQ(info);
	      
  // I don't need this
  info = PetscPrintf(PETSC_COMM_WORLD,"Program exit\n"); CHKERRQ(info);


  // Finalize  
  //  info = TaoRegisterDestroy(); CHKERRQ(info);
  info = TaoDestroy(tao);  CHKERRQ(info);
  info = TaoApplicationDestroy(taoapp); CHKERRQ(info);

  info = VecDestroy(X); CHKERRQ(info);
  info = VecDestroy(G); CHKERRQ(info);
  info = VecDestroy(X0); CHKERRQ(info);
  info = VecDestroy(L); CHKERRQ(info);
  info = VecDestroy(U); CHKERRQ(info);
  info = MatDestroy(H); CHKERRQ(info);
  

  info = TaoFinalize(); CHKERRQ(info);
  info = PetscFinalize(); CHKERRQ(info);


}


#if defined USING_SECOND_FUNCTIONS
//==========================================
#undef __FUNCT__
#define __FUNCT__ "funcgrad2"
int funcgrad2(TAO_SOLVER solver, TaoVec *XX, double *f, TaoVec *GG, void *ctx)
{
  double y1,y2,y3,A,B,C;
  Vec X,G;
  double *x, *g;
  int info;

  TaoFunctionBegin;

  // Get a handle to the current vector
  info = TaoVecGetPetscVec(XX,&X); CHKERRQ(info);
  info = VecGetArray(X,&x); CHKERRQ(info);
  
  // Get a handle to the current gradient
  info = TaoVecGetPetscVec(GG,&G); CHKERRQ(info);
  info = VecGetArray(G, &g); CHKERRQ(info);

  y1 = 1-x[1];
  y2 = 1-x[1]*x[1];
  y3 = 1-x[1]*x[1]*x[1];

  
  A = 1.5 - x[0]*y1;
  B = 2.25 - x[0]*y2;
  C = 2.625 - x[0]*y3;

  
  *f = A*A + B*B + C*C;

  g[0] = -2*(A*y1 + B*y2 + C*y3);
  g[1] = 2*x[0]*(A + 2*B*x[1] + 3*C*x[1]*x[1]);


  // restore the handles
  info = VecRestoreArray(X, &x); CHKERRQ(info);
  info = VecRestoreArray(G, &g); CHKERRQ(info);
  TaoFunctionReturn(0);
}

//==========================================
#undef __FUNCT__
#define __FUNCT__ "hessian2"
int hessian2(TAO_SOLVER solver, TaoVec *XX, TaoMat *HH, void *ctx)
{

  Mat H;
  PetscScalar v[4];
  double *xptr;
  Vec X;
  int idxm[2],idxn[2];
  int i;
  double f1,f2,f3,x,y;
  int info;

  TaoFunctionBegin;

  // Get a handle to the current vector
  info = TaoVecGetPetscVec(XX,&X); CHKERRQ(info);
  info = VecGetArray(X,&xptr); CHKERRQ(info);
  x = xptr[0];
  y = xptr[1];
  

  // Calculate the hessian, store in v


  f1 = 1.5 - x*(1-y);
  f2 = 2.25 - x*(1-y*y);
  f3 = 2.625 - x*(1-y*y*y);
  
  v[0] = 2*(y-1)*(y-1) + 2*(y*y-1)*(y*y-1) + 2*(y*y*y-1) * (y*y*y-1);
  v[3] = 2*x*x + 4*f2*x + 4*x*y*2*x*y + 12*f3*y*x + 6*y*y*x + 3*x*y*y;
  v[1] = v[2] = 2*f1 + 2*x*(y-1) + 4*f2*y + 4*y*x*(y*y-1) + 6*f3*y*y
    + 6*y*y*x*(y*y*y-1);
  


  /*
  v[0] = 6-4*x[1] - 2*x[1]*x[1] - 4*x[1]*x[1]*x[1] + 2*x[1]*x[1]*x[1]*x[1] + 2*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
  v[1] = v[2] = 3 - 4*x[0] - 4*x[0]*x[1] + 9*x[1] * 8*x[0]*x[1]*x[1]*x[1] + 15.75*x[1]*x[1] - 12*x[0]*x[1]*x[1] + 
    12*x[1]*x[1]*x[1]*x[1]*x[1];
  v[3] = -2*x[0]*x[0] + 9*x[0] +12*x[0]*x[0]*x[1]*x[1] + 31.5*x[0]*x[1] - 12*x[0]*x[0]*x[1] + 30*x[0]*x[0]*x[1]*x[1]*x[1]*x[1];
  */

  // Natural indices
  for (i=0;i<2;i++)
    idxm[i] = idxn[i] = i;

  info = TaoMatGetPetscMat(HH, &H); CHKERRQ(info);
  info = MatSetValues(H,2,idxm,2,idxn,v,INSERT_VALUES); CHKERRQ(info);

  info = VecRestoreArray(X,&xptr); CHKERRQ(info);


  info = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);


  // Display the hessian
  info = MatView(H,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(info);


  
  TaoFunctionReturn(0);
}

#else /* Using "first" functions */
//==========================================
#undef __FUNCT__
#define __FUNCT__ "funcgrad"
/* Function:  f(x,y) = (x-4)^2 + (y-3)^2 */
int funcgrad(TAO_SOLVER solver, Vec X, double *f, Vec G, void *ctx)
{

  double *x, *g;
  double ff;
  int i;
  int *ix;
  int info;
  AppCtx *user = (AppCtx *)ctx;

  TaoFunctionBegin;

  *f = user->c;

  // Retrieve the TaoVec object XX into the array x
  info =  VecGetArray(X, &x); CHKERRQ(info);



  //  Set the function value
  ff = (x[0]-4)*(x[0]-4) + (x[1]-3)*(x[1]-3);
  *f = ff;
  
  // Calculate the gradient
  info = VecGetArray(G, &g); CHKERRQ(info);

  g[0] = 2*(x[0]-4);
  g[1] = 2*(x[1]-3);


  info = VecRestoreArray(G, &g); CHKERRQ(info);
  info = VecRestoreArray(X,&x); CHKERRQ(info);

  TaoFunctionReturn(0);
}

//==========================================
#undef __FUNCT__
#define __FUNCT__ "hessian"
/* Stores the hessian of the function:  f(x,y) = (x-4)^2 + (y-3)^2  at the point XX */
int hessian(TAO_SOLVER solver, Vec X, Mat *H, Mat *Pre, MatStructure *flg, void *ctx)
{

  int im[2], in[2];
  int i;
  PetscScalar v[4];
  int info;


  TaoFunctionBegin;
  // Allocate space for the 2-d array and set the hessian values
  v[0] = v[3] = 2.0;
  v[1] = v[2] = 0.0;


  for (i=0;i<2;i++)
    im[i] = in[i] = i;

  // Now set the values
  info = MatSetValues(*H,2,im,2,in,v,INSERT_VALUES); CHKERRQ(info);

  info = MatAssemblyBegin(*H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(*H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);

  TaoFunctionReturn(0);
}

//==========================================
#undef __FUNCT__
#define __FUNCT__ "my_monitor"
int my_monitor(TAO_SOLVER solver, void *mctx)
{

  int iter;
  double f,gnorm,cnorm,xdiff;
  TaoTerminateReason reason;
  double *x;
  Vec X;
  TaoVec *XX;
  int rank;
  int info;


  TaoFunctionBegin;

  info = TaoGetSolution(solver,&XX); CHKERRQ(info);
  info = TaoVecGetPetscVec(XX,&X); CHKERRQ(info);
  info = VecGetArray(X,&x); CHKERRQ(info);

  
  info = PetscPrintf(PETSC_COMM_SELF,"x[0] = %f\nx[1] = %f\n",x[0],x[1]); CHKERRQ(info);


  info = TaoGetIterationData(solver, &iter, &f, &gnorm, &cnorm, &xdiff, &reason); CHKERRQ(info);
  
  info = PetscPrintf(PETSC_COMM_SELF,"iterations = %6d\tobj. value = %8.6g\tgnorm=%8.6g\n",iter,f,gnorm); CHKERRQ(info);

  if (reason)
  {
    info = PetscPrintf(PETSC_COMM_SELF,"Reason for termination is [%d]\n",reason); CHKERRQ(info);
  }

  info = VecRestoreArray(X,&x); CHKERRQ(info);
  
  TaoFunctionReturn(0);
  
}


#endif /* !USING_SECOND_FUNCTIONS */
