#include "petscda.h"
#include "tao.h"
#include "asl.h"
#include "getstub.h"

/* ampl variables -- stored in asl structure, but #defined to resemble globals: 
   int n_var      -- # of variables
   double LUv[]   -- lower bounds on variables
   double Uvx[]   -- upper bounds on variables
   double X0[]    -- initial vector (not defined if havex0=0)
   char *havex0   -- NULL if using default(0 vector)
   
*/
   
   
static char help[] = "";

/* Global variables */
ASL *asl;
int *indices; /* used to set PETSc vec's */

int FormFunctionGradient(TAO_APPLICATION,Vec,double *,Vec,void*);
int SetBounds(TAO_APPLICATION, Vec, Vec, void*);

double fatol,frtol,gatol,grtol,gttol,f_min;
int max_it,max_funcs,lmvm_lm,monitor;

//------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "check_problem_type"
/* Make sure it's a bound-constrained problem */
void check_problem_type()
{
  char error_message[1000]={0};
  int wrong = 0;
#define Assert(x) if(!(x)) {wrong=1; fprintf(stderr,#x" not satisfied!\n");}
  Assert(nbv == 0);  //number of binary variables
  Assert(niv == 0);  //# of other integer variables
  Assert(n_con == 0);//# of constraints
  Assert(n_obj == 1);//# of objectives
  Assert(nlo == 0 || nlo == 1);//# of nonlinear objectives
  Assert(nranges == 0);// number of ranged constraints
  Assert(nlc == 0);  //# of nonlinear general constraints
  Assert(nlnc == 0); //# of nonlinear network constraints
  Assert(nlvb == 0); //# of variables appearing nonlinearly in both con and obj
  Assert(nlvbi == 0);//# of integer variables appearing nonlinearly in both con and obj 
  Assert(nlvc == 0); //# of variables appearing nonlinearly in constraints
  Assert(nlvci == 0);//# of integer variables appearing nonlinearly in constraints
  Assert(nlvoi == 0);//# of integer variables appearing nonlinearly in objectives
  Assert(lnc == 0);  //# of linear network constraints
  if(wrong)
  {
    strcpy(error_message,"Error: Blmvm cannot handle this problem.  \n\
         Most likely, your model has constraints (n_con) other than variable bounds \n\
         or has integer (niv) or binary (nbv) variables.\n");
    write_sol(error_message,0,0,0);
    exit(0);	/* hack, or AMPL will ignore our output */
  }
}


//------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[])
{
  FILE *nl;
  Vec X;
  // initial, lower bound, upper bound, gradient vectors  
  double *x;
  TAO_SOLVER tao;
  TAO_APPLICATION taoapp;            /* TAO application context */
  TaoMethod method = "tao_blmvm";
  char lmvmstr[10];
  double f;  /* objective function value */
  int i,info,dim;
  PetscTruth flag;
  char *stub;
  TaoTerminateReason reason;

  /* AMPL-related variables */

  static keyword keywds[] = {      /* must be in alphabetical order */
    KW("tao_fatol",   D_val, &fatol,    "Absolute solution stopping tolerance."),
    KW("tao_fmin",    D_val, &f_min,     "Minimum function value stopping tolerance."),
    KW("tao_frtol",   D_val, &frtol,    "Relative solution stopping tolerance."),
    KW("tao_gatol",   D_val, &gatol,    "Absolute gradient norm stopping tolerance."),
    KW("tao_grtol",   D_val, &grtol,    "Relative gradient norm to function value stopping tolerance."),
    KW("tao_gttol",   D_val, &gttol,    "Gradient norm relative to initial times factor stopping tolerance."),
    KW("tao_lmvm_lm", I_val, &lmvm_lm,      "Vectors in Hessian approximation."),
    KW("tao_max_funcs", I_val, &max_funcs,      "Maximum function evaluations."),
    KW("tao_max_it",  I_val, &max_it,   "Maximum iterations."),
    KW("tao_monitor",  I_val, &monitor,   "Level of output."),
  };

  /*
    One side effect here is that options set by the user will
     be printed to stdout by AMPL.
  */
  static Option_Info Oinfo = {
    "blmvm", "BLMVM", "blmvm_options", keywds, nkeywds
  };

  info = PetscInitialize(&argc,&argv,(char *)0,help);
  info = TaoInitialize(&argc, &argv, (char *)0,help); CHKERRQ(info);

  /* Set default option values */
  fatol = 1.0e-12;
  frtol = 1.0e-20;
  gatol = 1.0e-6;
  grtol = 0.0;
  gttol = 0.0;
  max_it = 100000;
  max_funcs = 100000;
  f_min = -1e30;
  lmvm_lm = 5;
  monitor = 0;

  /* Read problem from ampl */
  asl = ASL_alloc(ASL_read_fg);  // partially separable function and gradient
  amplflag = 1;  
  solve_result_num = 500;  /* set failure value unless otherwise proven */
  stub = getstub(&argv, &Oinfo);
  if (!stub) {
    usage_ASL(&Oinfo, 1);
  }

  /* Get variable information */
  nl = jac0dim(stub,strlen(stub));

  /* Check to make sure problem is unconstrained, no integer variables */
  check_problem_type();
  
  /* Use AMPL function to read options */
  if (getopts(argv, &Oinfo))
    exit(2);

  /* Set the global vectors of asl */
  info = PetscMalloc(n_var*sizeof(real),(real**)&LUv); CHKERRQ(info);
  info = PetscMalloc(n_var*sizeof(real),(real**)&Uvx); CHKERRQ(info);
  info = PetscMalloc(n_var*sizeof(real),(real**)&X0); CHKERRQ(info);


  /* Set the initial vector to zero's */
  info = PetscMemzero(X0,n_var*sizeof(real));

  /* Read in the problem */
  fg_read(nl,0); 

  /* Ignore the -AMPL switch (AMPL adds this switch to the command line.
     If this program doesn't ask about the flag, PETSc will give a warning. */
  info = PetscOptionsHasName(0,"-AMPL",&flag); CHKERRQ(info);

  /* Create the vectors */
  info = VecCreateSeq(PETSC_COMM_SELF,n_var,&X); CHKERRQ(info);

  /* Create the Tao solver */
  info = TaoApplicationCreate(PETSC_COMM_WORLD,&taoapp); CHKERRQ(info);
  info = TaoCreate(PETSC_COMM_WORLD, method, &tao); CHKERRQ(info);

  info = TaoSetTolerances(tao,fatol,frtol,0,0); // The 0's are constraint tol's
  CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,gatol,grtol,gttol); CHKERRQ(info);
  info = TaoSetMaximumIterates(tao,max_it); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,max_funcs); CHKERRQ(info);
  info = TaoSetFunctionLowerBound(tao,f_min); CHKERRQ(info);
  snprintf(lmvmstr,9,"%d",lmvm_lm);
  info = PetscOptionsSetValue("-tao_lmm_vectors",lmvmstr); CHKERRQ(info);
  if (monitor)
  {
    info = TaoSetMonitor(tao,TaoDefaultMonitor,0); CHKERRQ(info);
  }


  /* Set objective function and gradient evaluation routines */
  info = TaoAppSetObjectiveAndGradientRoutine(taoapp,FormFunctionGradient,(void*)0); CHKERRQ(info);

  
  info = TaoAppSetVariableBoundsRoutine(taoapp,SetBounds,(void*)0); CHKERRQ(info);


  /* Set the initial vector (if necessary) */
  if (X0)
  {
  /* Set up the indices array.  indices[] is used to set the initial vector*/
    info = PetscMalloc(sizeof(int)*n_var,&indices); CHKERRQ(info);
    for (i=0;i<n_var;i++)
      indices[i] = i;

    info = VecSetValues(X,n_var,indices,X0,INSERT_VALUES); CHKERRQ(info);
    info = VecAssemblyBegin(X); CHKERRQ(info); 
    info = VecAssemblyEnd(X); CHKERRQ(info); 
    info = TaoAppSetInitialSolutionVec(taoapp, X); CHKERRQ(info);
    info = PetscFree(indices); CHKERRQ(info);
  }


  /* Solve the problem */


  info = TaoSolveApplication(taoapp, tao); CHKERRQ(info);

  info = TaoGetTerminationReason(tao,&reason); CHKERRQ(info);
  /* Output results */
  info = TaoView(tao); CHKERRQ(info);

  printf("\n-----------------------------------\n");

  /* return x to ampl (with message) */
  info = VecGetArray(X,&x); CHKERRQ(info); 
  if (reason>0) {
    solve_result_num = 100;
    write_sol("Solution found.",x,0,0); 
  }
  else
  {
    write_sol("Solution not found, returning last iteration.\n",x,0,0);
  }
  info = VecRestoreArray(X,&x); CHKERRQ(info); 


  /* Destroy the data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(taoapp);  CHKERRQ(info);
  info = VecDestroy(X); CHKERRQ(info);


  info = PetscFree(X0); CHKERRQ(info);
  info = PetscFree(LUv); CHKERRQ(info);
  info = PetscFree(Uvx); CHKERRQ(info);

  /* Finalize TAO */
  PetscFinalize();
  TaoFinalize();
  
  return 0;

}

//------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FormFunctionGradient"
int FormFunctionGradient(TAO_APPLICATION taoapp,Vec X,double *f,Vec G,void *)
{
  int info;
  double *x,*g;
  int i;

  /* Get handles to the vector arrays */
  info = VecGetArray(X,&x); CHKERRQ(info);
  info = VecGetArray(G,&g); CHKERRQ(info);
  

  /* Call the asl function and gradient evaluation routines */
  /* fint is #defined in asl.h */
  *f = objval(0,x,(fint*)&info); CHKERRQ(info);
  objgrd(0,x,g,(fint*)&info); CHKERRQ(info);

  /* return the handles */
  info = VecRestoreArray(X,&x); CHKERRQ(info);
  info = VecRestoreArray(G,&g); CHKERRQ(info);
  

  return 0;

}

int SetBounds(TAO_APPLICATION taoapp, Vec XL, Vec XU, void*)
{

  double *xl, *xu;
  int i, info;
  
  /* Set the upper and lower bounds */
  info = VecGetArray(XL,&xl); CHKERRQ(info);
  info = VecGetArray(XU,&xu); CHKERRQ(info);

  //  Lower and upper bounds stored in LUv and Uvx
  for (i=0;i<n_var;i++)
  {
    xl[i] = (LUv[i] > TAO_NINFINITY) ? LUv[i] : TAO_NINFINITY;
    xu[i] = (Uvx[i] < TAO_INFINITY) ? Uvx[i] : TAO_INFINITY;
  }

  // return the handles
  info = VecRestoreArray(XL,&xl); CHKERRQ(info);
  info = VecRestoreArray(XU,&xu); CHKERRQ(info);


  return 0;
}
