#include <iostream>
#include <string>
#include <stdlib.h>
#include "tao.h"
#include "Optimize_OptimizationModel.hxx"
#include "Solver_OptimizationSolver.hxx"
#include "Solver_ProjectState.hxx"
#include "sidl.hxx"
#include "sidl_DLL.hxx"

sidl::BaseClass LOADCLASS(std::string lib, std::string cls) ;

int main(int argc, char *argv[])
{
  // Create the sidl structures
  Optimize::OptimizationModel model= sidl::babel_cast<Optimize::OptimizationModel>(LOADCLASS("libRosenbrock-server-C++.so","Rosenbrock.RosenbrockModel"));
  Solver::ProjectState taoEnvironment = sidl::babel_cast<Solver::ProjectState>(LOADCLASS("libTaoapi-server-C++.so","TAO.Environment"));
  Solver::OptimizationSolver tao = sidl::babel_cast<Solver::OptimizationSolver>(LOADCLASS("libTaoapi-server-C++.so","TAO.Solver"));
  int info;

  // Initialize tao
  sidl::array<std::string> args = sidl::array<std::string>::create1d(argc);
  for (int i=0;i<argc;i++)
    args.set(i,argv[i]);

  info = taoEnvironment.Initialize(args); CHKERRQ(info);
  
  // create the solver
  info = tao.Create("tao_lmvm"); CHKERRQ(info);

  // read in any tao command line arguments
  info = tao.SetFromOptions(); CHKERRQ(info);

  // Solve the application
  info = tao.SolveApplication(model); CHKERRQ(info);

  // View the solver information
  info = tao.View(); CHKERRQ(info);
  

  // Cleanup
  info = tao.Destroy(); CHKERRQ(info);


  taoEnvironment.Finalize();

  return 0;
}


sidl::BaseClass LOADCLASS(std::string lib, std::string cls) 
{
  sidl::BaseClass base = 0;
  ::sidl::DLL mydll = ::sidl::Loader::findLibrary(cls,
						  "ior/impl",
						  ::sidl::Scope_SCLSCOPE, 
						  ::sidl::Resolve_SCLRESOLVE);
  if  (mydll._not_nil()) {
    base = mydll.createClass(cls);
  } else {
    std::cout << "Could not load library " << lib << std::endl;
    std::cout << "SIDL_DLL_PATH = " << getenv("SIDL_DLL_PATH") << std::endl;
  }

  return base;
}



