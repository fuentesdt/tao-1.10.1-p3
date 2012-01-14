#!/bin/env python
# PYTHONPATH should have: $(BABEL_ROOT)/lib/pythonX.X/site-packages,
#                          $(TAO_DIR)/lib/pythonX.X/site-packages
#SIDL_DLL_PATH should have: $(TAO_DIR)/lib/$(PETSC_ARCH)
import sys
import os
import sidl.Scope
import sidl.Resolve
import sidl.BaseClass
import sidl.BaseInterface
import sidl.Loader
import Solver.OptimizationSolver
import Solver.ProjectState
import Optimize.OptimizationModel
import Numeric


def LOADCLASS(classname):
  dll = sidl.Loader.findLibrary(classname,
                                "ior/impl",
                                sidl.Scope.SCLSCOPE,
                                sidl.Resolve.SCLRESOLVE)
  if dll:
    base = dll.createClass(classname)
  else:
    print ("Error loading class " + classname)
    print("\tSIDL_DLL_PATH = " + os.getenv("SIDL_DLL_PATH"))
    print("\tPYTHONPATH = " + os.getenv("PYTHONPATH"))
    print("\tLD_LIBRARY_PATH = " + os.getenv("LD_LIBRARY_PATH"))
    sys.exit(1)
    
  return base
  

class Rosenbrock_driver:

  def __init__(self):
    pass

  def go(self):
    
    model = Optimize.OptimizationModel.OptimizationModel(LOADCLASS("Rosenbrock.RosenbrockModel"))
    tao = Solver.OptimizationSolver.OptimizationSolver(LOADCLASS("TAO.Solver"))
    taoenv = Solver.ProjectState.ProjectState(LOADCLASS("TAO.Environment"))

    args = Numeric.zeros(len(sys.argv),'O')
    for i in range(len(sys.argv)):
      args[i] = sys.argv[i]

    taoenv.Initialize(args)
    tao.Create("tao_lmvm")
    tao.SetFromOptions()
    tao.SolveApplication(model)
    tao.View()
    tao.Destroy()
    taoenv.Finalize()

if __name__ == "__main__":

  driver = Rosenbrock_driver()
  driver.go()
  
