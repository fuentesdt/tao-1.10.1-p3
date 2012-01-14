! "$Id$";
!
!  Program usage: mpirun -np 1 rosenbrock1f [-help] [all TAO options] 
!
!  Description:  This example demonstrates use of the TAO package to solve an
!  unconstrained minimization problem on a single processor.  We minimize the 
!  extended Rosenbrock function: 
!       sum_{i=0}^{n/2-1} ( alpha*(x_{2i+1}-x_{2i}^2)^2 + (1-x_{2i})^2 )
!
!  The C version of this code is rosenbrock1.c
!
!/*T
!  Concepts: TAO - Solving an unconstrained minimization problem
!  Routines: TaoInitialize(); TaoFinalize(); TaoSetOptions();
!  Routines: TaoApplicationCreate();
!  Routines: TaoCreate(); TaoAppSetObjectiveAndGradientRoutine(); 
!  Routines: TaoAppSetHessianMat(); TaoAppSetHessianRoutine();
!  Routines: TaoAppSetInitialSolutionVec(); 
!  Routines: TaoSolveApplication(); TaoDestroy(); TaoAppDestroy();
!  Routines: TaoView(); TaoGetTerminationReason();                              
!  Processors: 1
!T*/ 
!
! ---------------------------------------------------------------------- 
!
      implicit none

#include "rosenbrock1f.h"

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                   Variable declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!  See additional variable declarations in the file rosenbrock1f.h

      integer          info    ! used to check for functions returning nonzeros
      Vec              x       ! solution vector
      Mat              H       ! hessian matrix
      TAO_SOLVER       tao     ! TAO_SOVER context
      TAO_APPLICATION  taoapp   ! TAO application context
      PetscTruth       flg     
      PetscInt         i2,i1
      integer          size,rank    ! number of processes running
      PetscScalar      zero
      TaoTerminateReason reason

      


!  Note: Any user-defined Fortran routines (such as FormGradient)
!  MUST be declared as external.

      external FormFunctionGradient,FormHessian

      zero = 0.0d0
      i2 = 2
      i1 = 1

!  Initialize TAO and PETSc
      call PetscInitialize(PETSC_NULL_CHARACTER,info)
      call TaoInitialize(PETSC_NULL_CHARACTER,info)

      call MPI_Comm_size(PETSC_COMM_WORLD,size,info)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,info)
      if (size .ne. 1) then
         if (rank .eq. 0) then
            write(6,*) 'This is a uniprocessor example only!'
         endif
         SETERRQ(1,' ',info)
      endif

!  Initialize problem parameters
      n     = 2
      alpha = 99.0d0



! Check for command line arguments to override defaults        
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,           &
     &                        info)
      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-alpha',            &
     &                           alpha,flg,info)

!  Allocate vectors for the solution and gradient
      call VecCreateSeq(PETSC_COMM_SELF,n,x,info)

!  Allocate storage space for Hessian; 
      call MatCreateSeqBAIJ(PETSC_COMM_SELF,i2,n,n,i1,                   &
     &     PETSC_NULL_INTEGER, H,info)

      call MatSetOption(H,MAT_SYMMETRIC,PETSC_TRUE,info)


!  The TAO code begins here 

!  Create TAO solver
      call TaoCreate(PETSC_COMM_SELF,TAO_NULL_CHARACTER,tao,info)
      call TaoApplicationCreate(PETSC_COMM_SELF,taoapp,info)
!  Set routines for function, gradient, and hessian evaluation 

!     TaoAppSetObjectiveAndGradientRoutine is shortened to 31 chars to comply with some compilers
      call TaoAppSetObjectiveAndGradientRo(taoapp,                      &
     &     FormFunctionGradient,TAO_NULL_OBJECT,info)
      call TaoAppSetHessianMat(taoapp,H,H,info)
      call TaoAppSetHessianRoutine(taoapp,FormHessian,TAO_NULL_OBJECT,   &
     &     info)

!  Optional: Set initial guess
      call VecSet(x, zero, info)
      call TaoAppSetInitialSolutionVec(taoapp, x, info)


!  Check for TAO command line options
      call TaoSetTolerances(tao,1.0d-3,1.0d-3,1.0d-3,1.0d-3,info)
      call TaoSetOptions(taoapp,tao,info)

!  SOLVE THE APPLICATION
      call TaoSolveApplication(taoapp,tao,info)

      call TaoGetTerminationReason(tao, reason, info)
      if (reason .le. 0) then
         print *,'Try a different TAO method, adjust some parameters,'
         print *,'or check the function evaluation routines.'
      endif

!  TaoView() prints info about the TAO solver; the option
!      -tao_view 
!  can alternatively be used to activate this at runtime.
!     call TaoView(tao,info)
      

!  Free TAO data structures
      call TaoDestroy(tao,info)
      call TaoAppDestroy(taoapp,info)

!  Free PETSc data structures 
      call VecDestroy(x,info)
      call MatDestroy(H,info)

!  Finalize TAO 
      call TaoFinalize(info)
      call PetscFinalize(info)

      end


! --------------------------------------------------------------------
!  FormFunctionGradient - Evaluates the function f(X) and gradient G(X)
!
!  Input Parameters:
!  tao - the TAO_SOLVER context
!  X   - input vector
!  dummy - not used
!
!  Output Parameters:
!  G - vector containing the newly evaluated gradient
!  f - function value
      
      subroutine FormFunctionGradient(taoapp, X, f, G, dummy, info)
      implicit none

! n,alpha defined in rosenbrock1f.h
#include "rosenbrock1f.h"

      TAO_APPLICATION  taoapp
      Vec              X,G
      PetscScalar      f
      integer          info
      PetscInt         dummy


      PetscScalar      ff,t1,t2
      PetscInt         i,nn

! PETSc's VecGetArray acts differently in Fortran than it does in C.
! Calling VecGetArray((Vec) X, (PetscScalar) x_array(0:1), (PetscOffset) x_index, info)
! will return an array of doubles referenced by x_array offset by x_index.
!  i.e.,  to reference the kth element of X, use x_array(k + x_index).
! Notice that by declaring the arrays with range (0:1), we are using the C 0-indexing practice.
      PetscScalar      g_v(0:1),x_v(0:1)
      PetscOffset      g_i,x_i

      info = 0
      nn = n/2
      ff = 0

!     Get pointers to vector data
      call VecGetArray(X,x_v,x_i,info)
      call VecGetArray(G,g_v,g_i,info)


!     Compute G(X)
      do i=0,nn-1
         t1 = x_v(x_i+2*i+1) - x_v(x_i+2*i)*x_v(x_i+2*i)
         t2 = 1.0 - x_v(x_i + 2*i)
         ff = ff + alpha*t1*t1 + t2*t2
         g_v(g_i + 2*i) = -4*alpha*t1*x_v(x_i + 2*i) - 2.0*t2
         g_v(g_i + 2*i + 1) = 2.0*alpha*t1
      enddo

!     Restore vectors
      call VecRestoreArray(X,x_v,x_i,info)
      call VecRestoreArray(G,g_v,g_i,info)

      f = ff
      call PetscLogFlops(nn*15,info)

      return
      end

!  
! ---------------------------------------------------------------------
!
!  FormHessian - Evaluates Hessian matrix.
!
!  Input Parameters:
!  tao     - the TAO_SOLVER context
!  X       - input vector
!  dummy   - optional user-defined context, as set by SNESSetHessian()
!            (not used here)
!
!  Output Parameters:
!  H      - Hessian matrix
!  PrecH  - optionally different preconditioning matrix (not used here)
!  flag   - flag indicating matrix structure
!  info   - error code
!
!  Note: Providing the Hessian may not be necessary.  Only some solvers
!  require this matrix.

      subroutine FormHessian(taoapp,X,H,PrecH,flag,dummy,info)
      implicit none

#include "rosenbrock1f.h"

!  Input/output variables:
      TAO_APPLICATION  taoapp
      Vec              X
      Mat              H, PrecH
      MatStructure     flag
      integer          info
      PetscInt         dummy
      
      PetscScalar      v(0:1,0:1)
      PetscTruth assembled

! PETSc's VecGetArray acts differently in Fortran than it does in C.
! Calling VecGetArray((Vec) X, (PetscScalar) x_array(0:1), (PetscOffset) x_index, info)
! will return an array of doubles referenced by x_array offset by x_index.
!  i.e.,  to reference the kth element of X, use x_array(k + x_index).
! Notice that by declaring the arrays with range (0:1), we are using the C 0-indexing practice.
      PetscScalar      x_v(0:1)
      PetscOffset      x_i
      PetscInt         i,nn,ind(0:1),i2


      info = 0
      nn= n/2
      i2 = 2

!  Zero existing matrix entries
      call MatAssembled(H,assembled,info)
      if (assembled .eqv. PETSC_TRUE) call MatZeroEntries(H,info)

!  Get a pointer to vector data

      call VecGetArray(X,x_v,x_i,info)

!  Compute Hessian entries

      do i=0,nn-1
         v(1,1) = 2.0*alpha
         v(0,0) = -4.0*alpha*(x_v(x_i+2*i+1) -                          &
     &                3*x_v(x_i+2*i)*x_v(x_i+2*i))+2
         v(1,0) = -4.0*alpha*x_v(x_i+2*i)
         v(0,1) = v(1,0)
         ind(0) = 2*i
         ind(1) = 2*i + 1
         call MatSetValues(H,i2,ind,i2,ind,v,INSERT_VALUES,info)
      enddo

!  Restore vector

      call VecRestoreArray(X,x_v,x_i,info)

!  Assemble matrix

      call MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY,info)
      call MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY,info)

      call PetscLogFlops(nn*9,info)

!  Set flag to indicate that the Hessian matrix retains an identical
!  nonzero structure throughout all nonlinear iterations (although the
!  values of the entries change). Thus, we can save some work in setting
!  up the preconditioner (e.g., no need to redo symbolic factorization for
!  ICC preconditioners).
!   - If the nonzero structure of the matrix is different during
!     successive linear solves, then the flag DIFFERENT_NONZERO_PATTERN
!     must be used instead.  If you are unsure whether the matrix
!     structure has changed or not, use the flag DIFFERENT_NONZERO_PATTERN.
!   - Caution:  If you specify SAME_NONZERO_PATTERN, the software
!     believes your assertion and does not check the structure
!     of the matrix.  If you erroneously claim that the structure
!     is the same when it actually is not, the new preconditioner
!     will not function correctly.  Thus, use this optimization
!     feature with caution!

      flag = SAME_NONZERO_PATTERN

      return
      end




