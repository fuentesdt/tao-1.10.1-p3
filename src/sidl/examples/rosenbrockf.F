      program testdriver
      implicit none
     
      integer*8 model_base, model
      integer*8 tao_base, tao
      integer*8 taoenv_base, taoenv
      integer*8 args
      integer info,i
      character*80 arg


c *** Load the classes from SIDL_DLL_PATH
      call sidl_DLL_createClass_f('libRosenbrock-server-C++.so',        &
     &                            'Rosenbrock.RosenbrockModel',         &
     &                             model_base)
      if (model_base.eq.0) then
         print *,"unable to load class Rosenbrock.RosenbrockModel"
      endif

      call sidl_DLL_createClass_f('libTaoapi-server-C++.so',            &
     &                            'TAO.Solver',tao_base)
      if (tao_base.eq.0) then
         print *,"unable to load class TAO.Solver"
      endif

      call sidl_DLL_createClass_f('libTaoapi-server-C++.so',            &
     &                            'TAO.Environment',taoenv_base)
      if (taoenv_base .eq. 0) then
         print *,"unable to load class TAO.Environment"
      endif

c *** Cast to usable interfaces
      call Optimize_OptimizationModel__cast_f(model_base, model)
      call Solver_OptimizationSolver__cast_f(tao_base, tao)
      call Solver_ProjectState__cast_f(taoenv_base, taoenv)


c *** Get command line arguments
      call sidl_string__array_create1d_f(args,iargc()+1)
      call sidl_string__array_set1_f(args,0,'rosenbrockf')
      do i=1,iargc()
         call getarg(i,arg)
         call sidl_string__array_set1_f(args,i,arg)
      enddo
c *** Initializes PETSc and TAO
      call Solver_ProjectState_Initialize_f(taoenv,args,info)

c ***  Create the solver 
      call Solver_OptimizationSolver_Create_f(tao,"tao_lmvm",info)

c *** Read in any tao command line arguments
      call Solver_OptimizationSolver_SetFromOptions_f(tao,info)

c *** Solve the application 
      call Solver_OptimizationSolver_SolveApplication_f(tao,model,info)

c *** View the solver information 
      call Solver_OptimizationSolver_View_f(tao,info)


      

c *** delete sidl structures
      call Optimize_OptimizationModel_deleteRef_f(model)
      call Solver_OptimizationSolver_deleteRef_f(tao)
      call Solver_ProjectState_Finalize_f(taoenv,info)

      end

      

      

      
      
