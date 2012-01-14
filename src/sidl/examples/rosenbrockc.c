//#include "tao.h"
#include <stdio.h>
#include <stdlib.h>
#include "Optimize_OptimizationModel.h"
#include "Solver_OptimizationSolver.h"
#include "Solver_ProjectState.h"
#include "sidl_BaseClass.h"
#include "sidl.h"

void CHKERR(int info, char*);
sidl_BaseClass LOADCLASS(char *, char *);

int main(int argc, char *argv[]) 
{
  /* Create the sidl structures */
  sidl_BaseClass model_base = LOADCLASS("libRosenbrock-server-C++.so","Rosenbrock.RosenbrockModel");
  Optimize_OptimizationModel model = Optimize_OptimizationModel__cast(model_base);

  sidl_BaseClass tao_base = LOADCLASS("libTaoapi-server-C++.so","TAO.Solver");
  Solver_OptimizationSolver tao = Solver_OptimizationSolver__cast(tao_base);

  sidl_BaseClass taoenvironment_base = LOADCLASS("libTaoapi-server-C++.so","TAO.Environment");
  Solver_ProjectState taoenvironment = Solver_ProjectState__cast(taoenvironment_base);


  int info,i;

  /* Get argument list in sidl_string__array form */
  struct sidl_string__array *args = sidl_string__array_create1d(argc);
  for (i=0;i<argc;i++) 
    sidl_string__array_set1(args, i, argv[i]);


  /* Initializes PETSc and TAO */
  info = Solver_ProjectState_Initialize(taoenvironment,args); CHKERR(info, "initialize");

  /* Create the solver */
  info = Solver_OptimizationSolver_Create(tao,"tao_lmvm");   CHKERR(info,"solver_create");

  /* Read in any tao command line arguments */
  info = Solver_OptimizationSolver_SetFromOptions(tao);   CHKERR(info,"setfromoptions");

  /* Solve the application */
  info = Solver_OptimizationSolver_SolveApplication(tao,model);   CHKERR(info,"solveapplication");

  /* View the solver information */
  info = Solver_OptimizationSolver_View(tao);   CHKERR(info,"view");


  /* Clean up */
  info = Solver_OptimizationSolver_Destroy(tao);   CHKERR(info,"solver_destroy");
  info = Solver_ProjectState_Finalize(taoenvironment); CHKERR(info,"env_destroy");



  /* Delete sidl structure */
  Optimize_OptimizationModel_deleteRef(model);
  Solver_OptimizationSolver_deleteRef(tao);
  sidl_string__array_deleteRef(args);
  return 0;
}


void CHKERR(int info, char *str)
{
  if (info) {
    printf("ERROR: %d: %s\n",info,str);
    exit(info);
  }
  return;
}

sidl_BaseClass LOADCLASS(char *lib, char *cls) 
{
  sidl_BaseClass base=0;
  sidl_DLL dll = sidl_Loader_findLibrary(lib,
					 "ior/impl",
					 sidl_Scope_SCLSCOPE,
					 sidl_Resolve_SCLRESOLVE);
  if (dll) {
    base = sidl_DLL_createClass(dll,cls);
  } else {
    fprintf(stderr,"Could not load library %s\n",lib);
    fprintf(stderr,"SIDL_DLL_PATH = %s\n",getenv("SIDL_DLL_PATH"));
    exit(1);
  }
  return base;
}
