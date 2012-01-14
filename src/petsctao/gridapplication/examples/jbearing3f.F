! Program usage: mpirun -np <proc> jbearing [-help] [all TAO options]
!

!/*T
!   Concepts: TAO - Solving a bounded minimization problem
!   Routines: TaoInitialize(); TaoFinalize();
!   Routines: TaoCreate(); TaoDestroy();
!   Routines: DAApplicationCreate(); DAApplicationDestroy();
!   Routines: DAAppSetVariableBoundsRoutine();
!   Routines: DAAppSetElementObjectiveAndGradientRoutine();
!   Routines: DAAppSetElementHessianRoutine();
!   Routines: DAAppSetObjectiveAndGradientRoutine();
!   Routines: DAAppSetHessianRoutine();
!   Routines: TaoSetOptions();
!   Routines: TaoGetIterationData(); TaoDAAppSolve();
!   Routines: DAAppSetBeforeMonitor(); DAAppSetAfterMonitor();
!   Routines: TaoView();
!   Routines: DAAppGetSolution();
!   Routines: DAAppGetInterpolationMatrix();
!   Processors: n
!T*/
      implicit none

#include "jbearing3f.h"

      integer           info         ! Error check
      PetscInt          iter         ! # of iterations solve took
      PetscInt          nlevels      ! multigrid levels
      PetscInt          Nx, Ny
      PetscScalar        ff, gnorm
      DA                 DAarray(0:19)
      Vec                x            ! solution vector
      PetscTruth         flg
      TAO_SOLVER         tao
      TAO_APPLICATION    JBearApp
      TaoTerminateReason reason
      PetscInt          i

      external DASetBounds, MyGridMonitorBefore, MyGridMonitorAfter
      external JBearLocalHessian, JBearLocalFunctionGradient
      external WholeJBearHessian, WholeJBearFunctionGradient

! Inititialize TAO
      call PetscInitialize(PETSC_NULL_CHARACTER, info)
      call TaoInitialize(PETSC_NULL_CHARACTER, info)

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, info)

      ecc = 0.1d0
      b = 10.0d0
      mx = 11
      my = 11
      nlevels = 5
      
      call PetscOptionsGetReal(TAO_NULL_CHARACTER,"-ecc",ecc,flg,info)
      call PetscOptionsGetReal(TAO_NULL_CHARACTER,"-b",b,flg,info)
      call PetscOptionsGetInt(TAO_NULL_CHARACTER,"-mx",mx,flg,info)
      call PetscOptionsGetInt(TAO_NULL_CHARACTER,"-my",my,flg,info)
      call PetscOptionsGetInt(TAO_NULL_CHARACTER,"-nlevels",nlevels,flg, &
     &     info)

! Let PETSc determine the vector distribution
      Nx = PETSC_DECIDE
      Ny = PETSC_DECIDE

! Create distributed array to manage parallel grid and vector
      call DACreate2d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,  &
     &     mx, my, Nx, Ny, 1, 1, PETSC_NULL_INTEGER,                    &
     &     PETSC_NULL_INTEGER, DAarray(0), info)
      
      do i=1,nlevels-1
         call DARefine(DAarray(i-1), PETSC_COMM_WORLD, DAarray(i), info)
      enddo


! Create TAO solver and set desired solution method
      call TaoCreate(PETSC_COMM_WORLD, "tao_gpcg", tao, info)
      call TaoApplicationCreate(PETSC_COMM_WORLD, JBearApp, info)
      call TaoAppSetDAApp(JBearApp, DAarray, nlevels, info)

! Sets routine bounds evaluation
      call DAAppSetVariableBoundsRoutine(JBearApp, DASetBounds,          &
     &     TAO_NULL_OBJECT, info)

      call PetscOptionsHasName(TAO_NULL_CHARACTER,"-byelement",flg,info)
      if (flg .eqv. TAO_TRUE) then
!        Set routines for function and gradient evaluation, element by element

!     DaAppSetElementObjectiveAndGradientRoutine is shortened to comply with some compilers
         call DAAppSetElementObjectiveAndGrad(JBearApp,                  &
     &        JBearLocalFunctionGradient,63,TAO_NULL_OBJECT,info)
         call DAAppSetElementHessianRoutine(JBearApp, JBearLocalHessian, &
     &        16, TAO_NULL_OBJECT,info)
      else

!     DaAppSetObjectiveAndGradientRoutine is shortened to comply with some compilers
         call DAAppSetObjectiveAndGradientRou(JBearApp,                  &
     &        WholeJBearFunctionGradient, TAO_NULL_OBJECT,info)
         call DAAppSetHessianRoutine(JBearApp, WholeJBearHessian,        &
     &        TAO_NULL_OBJECT, info)
      endif

      call DAAppSetBeforeMonitor(JBearApp, MyGridMonitorBefore,          &
     &     TAO_NULL_OBJECT, info)


! Check for -tao_monitor flag
      call PetscOptionsHasName(TAO_NULL_CHARACTER,"-tao_monitor",flg,    &
     &     info)

      if (flg .eqv. TAO_TRUE) then
         call DAAppSetAfterMonitor(JBearApp, MyGridMonitorAfter,         &
     &        TAO_NULL_OBJECT, info)
      endif


! Check for any tao command line options
      call TaoSetOptions(JBearApp, tao, info)

! SOLVE THE APPLICATION
      call TaoDAAppSolve(JBearApp, tao, info)

! Get information on terminaton
      call TaoGetSolutionStatus(tao, iter, ff, gnorm, TAO_NULL_SCALAR,   &
     &     TAO_NULL_SCALAR, reason, info)
      if (reason .le. 0 .and. rank .eq. 0) then
         print *,'Try a different TAO method, adjust some parameters,'
         print *,'or check the function evaluation routines.'
         print 110,iter,ff,gnorm
 110     format('Iterations:',I10,' Function Value:',f10.7,' Residual:',   &
     &        F10.7)
      endif


! Check for -view_sol flag
      call PetscOptionsHasName(TAO_NULL_CHARACTER,"-view_sol", flg,info)
      if (flg .eqv. TAO_TRUE) then
         call DAAppGetSolution(JBearApp, nlevels-1, x, info)
         call VecView(x, PETSC_VIEWER_STDOUT_WORLD,info)
      endif

! Free TAO data structures
      call TaoDestroy(tao, info)
      call TaoAppDestroy(JBearApp, info)

! Free PETSc data structures
      do i=0,nlevels-1
         call DADestroy(DAarray(i), info)
      enddo

      call TaoFinalize(info)
      call PetscFinalize(info)

      end

! Monitors to display info of iterations at every grid level

      subroutine MyGridMonitorBefore(myapp, da, level, dummy, info)
      implicit none
#include "jbearing3f.h"
      
      TAO_APPLICATION        myapp
      DA                     da
      PetscInt               level, dummy
      integer                info

      PetscScalar            t
      PetscInt               i

      call DAGetInfo(da, PETSC_NULL_INTEGER, mx, my, PETSC_NULL_INTEGER, &
     &     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,   &
     &     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,   &
     &     PETSC_NULL_INTEGER, info)

      hx = (2.0d0 * 3.14159265358979d0) / (mx - 1)
      hy = (2.0d0 * b) / (my - 1)
      area = 0.5d0 * hx * hy

      do i=0,mx-1
         t = 1.0d0 + ecc * cos(i*hx)
         wq(i) = t*t*t
         wl(i) = ecc * sin(i*hx)
      enddo

      call PetscLogFlops(8 + 7*mx,info)

      if (rank .eq. 0) print 100, level, mx, my
 100  format('Grid:',I6,' mx: ',I6, ' my: ', I6)

      info = 0
      end

      subroutine MyGridMonitorAfter(myapp, da, level, dummy, info)
      implicit none
#include "jbearing3f.h"

      TAO_APPLICATION myapp
      DA              da
      integer         info
      PetscInt        level, dummy

      PetscInt        n
      PetscScalar     dd
      Vec             X, XCoarse, XX
      Mat             Interpolate

      dd = -1.0d0

      call DAGetInfo(da, PETSC_NULL_INTEGER, mx, my, PETSC_NULL_INTEGER, &
     &     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,   &
     &     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,   &
     &     PETSC_NULL_INTEGER, info)


      if (level .gt. 0) then
         call DAAppGetInterpolationMatrix(myapp, level, Interpolate,     &
     &        TAO_NULL_OBJECT, info)
         call DAAppGetSolution(myapp, level, X, info)
         call DAAppGetSolution(myapp, level-1, XCoarse, info)
         call VecDuplicate (X,XX,info)
         call MatMult(Interpolate, XCoarse, XX, info)
         call VecAXPY(XX, dd, X, info)
         call VecNorm(XX, NORM_INFINITY, dd, info)
         if (rank .eq. 0) print 120,dd
 120     format('Maximum Interpolation Error: ',F10.7)
         call VecNorm(XX, NORM_1, dd, info)
         call VecGetSize(XX, n, info)
         if (rank .eq. 0) print 130,dd/n
 130     format('Average Interpolation Error: ',F10.7)
         call VecDestroy(XX, info)
      endif

      info = 0
      end


      subroutine DASetBounds(daapp, da, XL, XU, dummy, info)
      implicit none
#include "jbearing3f.h"

      TAO_APPLICATION     daapp
      DA                  da
      Vec                 XL, XU
      PetscInt            dummy
      integer             info

      PetscScalar         zero, large

      PetscScalar         xu_a(0:1)
      PetscOffset         xu_i

      PetscInt            xs, ys, xm, ym  ! corner information
      PetscInt            gxs,gys,gxm,gym ! ghost corners
      PetscInt            i,j,k

      Vec                 localXU
      
      PetscScalar         tmy, tmx

      zero = 0.0d0
      large = 1.0d20
      tmy = my
      tmx = mx

      call VecSet(XL, zero, info)
      call DAGetCorners(da, xs, ys, TAO_NULL_INTEGER, xm, ym,           &
     &     TAO_NULL_INTEGER, info)
      call DAGetGhostCorners(da, gxs, gys, TAO_NULL_INTEGER, gxm, gym,  &
     &     TAO_NULL_INTEGER, info)
      
      call DAGetLocalVector(da, localXU, info)
      call VecSet(localXU, large, info)
      call VecGetArray(localXU,xu_a,xu_i,info)
! xl[i,j] (in C)  = xl[i+(xm-1)*j] (in C) = xl_a(xl_i + (my-1)*i + j) in Fortran    

      do j = ys, ys+ym-1
         do i = xs, xs+xm-1
            k = i-gxs + (j-gys)*gxm
            if ((i .eq. 0) .or. (j .eq. 0) .or.                         &
     &           (i .eq. mx-1) .or. (j .eq. my-1)) then      
               xu_a(xu_i + k) = 0.0d0
            endif
         enddo
      enddo

      call VecRestoreArray(localXU,xu_a, xu_i, info)
      call DALocalToGlobal(da, localXU, INSERT_VALUES, XU, info)
      call DARestoreLocalVector(da,localXU,info)

      info = 0
      end


      subroutine JBearLocalFunctionGradient(coor, x, f, g, dummy, info)
      implicit none
#include "jbearing3f.h"

      PetscInt    coor(0:1)
      PetscScalar x(0:3),f,g(0:3)
      PetscInt    dummy
      integer     info

      PetscScalar avgWq, sqGrad, avgWV, fl, fu
      PetscScalar aread3, dvdx, dvdy
      PetscInt    i

      aread3 = area / 3.0d0
      i = coor(0)
      
! lower triangle contribution
      dvdx = (x(0) - x(1)) / hx
      dvdy = (x(0) - x(2)) / hy
      sqGrad = dvdx * dvdx + dvdy * dvdy
      avgWq = (2.0d0 * wq(i) + wq(i+1)) / 6.0d0
      avgWV = (wl(i) * x(0) + wl(i+1)*x(1) + wl(i)*x(2)) / 3.0d0
      fl = avgWq * sqGrad - avgWV

      dvdx = dvdx * hy * avgWq
      dvdy = dvdy * hx * avgWq
      
      g(0) = (dvdx + dvdy) - wl(i) * aread3
      g(1) = (-dvdx) - wl(i+1) * aread3
      g(2) = (-dvdy) - wl(i) * aread3

! upper triangle contribution
      dvdx = (x(3) - x(2)) / hx
      dvdy = (x(3) - x(1)) / hy
      sqGrad = dvdx * dvdx + dvdy * dvdy
      avgWq = (2.0d0 * wq(i+1) + wq(i)) / 6.0d0
      avgWV = (wl(i+1)*x(1) + wl(i)*x(2) + wl(i+1)*x(3)) / 3.0d0
      fu = avgWq * sqGrad - avgWV

      dvdx = dvdx * hy * avgWq
      dvdy = dvdy * hx * avgWq
      g(1) = g(1) - dvdy - wl(i+1)*aread3
      g(2) = g(2) - dvdx - wl(i)*aread3
      g(3) = (dvdx + dvdy) - wl(i+1) * aread3

      f = area * (fl + fu)
      info = 0
      end
      

      subroutine JBearLocalHessian(coor, x, H, dummy, info)
      implicit none
#include "jbearing3f.h"

      integer      coor(0:1)
      PetscScalar  x(0:3)
      PetscScalar  H(0:3,0:3)
      PetscInt     dummy
      integer      info

      PetscScalar  wql, wqu
      PetscScalar  dydx, dxdy, wqldydx
      PetscScalar  wqldxdy, wqudydx, wqudxdy
      PetscInt    i

      i = coor(0)

      dxdy = hx / hy
      dydx = hy / hx
      wql = (2.0d0 * wq(i) + wq(i+1)) / 6.0d0
      wqu = (wq(i) + 2.0d0 * wq(i+1)) / 6.0d0

      wqldydx = wql * dydx
      wqldxdy = wql * dxdy

      wqudydx = wqu * dydx
      wqudxdy = wqu * dxdy

! Hessian contribution at 0,0
      H(0,0) = wqldxdy + wqldydx
      H(1,0) = -wqldydx
      H(0,1) = -wqldydx
      H(0,2) = -wqldxdy
      H(2,0) = -wqldxdy
      H(3,0) = 0.0d0
      h(0,3) = 0.0d0
      

! Hessian contribution at 1,0
      H(1,1) = wqldydx + wqudxdy
      H(1,2) = 0.0d0
      H(2,1) = 0.0d0
      H(1,3) = -wqudxdy
      H(3,1) = -wqudxdy

! Hessian contribution at 0,2
      H(2,2) = wqldxdy + wqudydx
      H(2,3) = -wqudydx
      H(3,2) = -wqudydx

! Hessian contribution at 1,1
      H(3,3) = wqudydx + wqudxdy

      end


      subroutine WholeJBearFunctionGradient(daapp,da,X,f,G,dummy,info)
      implicit none
#include "jbearing3f.h"

      TAO_APPLICATION  daapp
      DA               da
      Vec              X
      PetscScalar      f
      Vec              G
      PetscInt         dummy
      integer          info

      Vec              localX
      PetscInt         i,j
      PetscInt         xs, ys, xm, ym
      PetscInt         xe, ye 
      PetscInt         gxs, gys, gxm, gym
      PetscInt         k, yep, xep, xsm, ysm
      PetscScalar      floc, fusum, flsum
      PetscScalar      zero
      PetscScalar      avgWq, sqGrad, avgWV
      PetscScalar      fl, fu, aread3
      PetscScalar      dvdx, dvdy
      PetscScalar      lx_a(0:1)            !local vector
      PetscOffset      lx_i
      PetscScalar      v, vt, vr, vl, vb    !triangle vertices
      PetscScalar      val1, val2, val3

      aread3 = area / 3.0d0
      zero = 0.0d0
      flsum = 0.0d0
      fusum = 0.0d0

      call VecSet(G, zero, info)
      
      call DAGetLocalVector(da, localX, info)

      call DAGlobalToLocalBegin(da, X, INSERT_VALUES, localX, info)
      call DAGlobalToLocalEnd(da, X, INSERT_VALUES, localX, info)

      call VecGetArray(localX, lx_a, lx_i, info)

      call DAGetCorners(da, xs, ys, TAO_NULL_INTEGER, xm, ym,           &
     &     TAO_NULL_INTEGER, info)
      call DAGetGhostCorners(da, gxs, gys, TAO_NULL_INTEGER, gxm, gym,  &
     &     TAO_NULL_INTEGER, info)

      xe = xs + xm 
      ye = ys + ym 


! Check to see if we are at an edge of the DA (if so, artificially set border)
      if (xs .eq. 0) then
         xsm = xs - 1
      else
         xsm = xs
      endif

      if (ys .eq. 0) then
         ysm = ys - 1
      else
         ysm = ys
      endif

      if (xe .eq. mx) then
         xep = xe + 1
      else
         xep = xe 
      endif
      
      if (ye .eq. my) then
         yep = ye + 1
      else
         yep = ye
      endif

      do j = ys, gys+gym-2       
         do i = xs, gxs + gxm - 2
            k = (j-gys)*gxm + i-gxs            ! index in local vector for this i,j

            ! lower triangle contribution
            v = zero
            vr = zero
            vt = zero
            fl = zero

            v = lx_a(lx_i + k)
            vr = lx_a(lx_i + k + 1)
            vt = lx_a(lx_i + k + gxm)
            dvdx = (v-vr)/hx
            dvdy = (v-vt)/hy

            sqGrad = dvdx * dvdx + dvdy * dvdy
            avgWq = (2.0d0 * wq(i) + wq(i+1)) / 6.0d0
            avgWV = (wl(i) * v +  wl(i+1)*vr + wl(i)*vt)/3.0d0
            fl = avgWq * sqGrad - avgWV

            dvdx = dvdx * hy * avgWq
            dvdy = dvdy * hx * avgWq

            val1 = dvdx + dvdy - wl(i) * aread3
            call VecSetValuesLocal(G,1,k,val1,ADD_VALUES,info)

            val2 = -dvdx - wl(i+1) * aread3
            call VecSetValuesLocal(G,1,k+1,val2,ADD_VALUES,info)

            val3 = -dvdy - wl(i)*aread3
            call VecSetValuesLocal(G,1,k+gxm,val3,ADD_VALUES,info)

            flsum = flsum + fl

            vb = zero
            vl = zero
            v = zero
            fu = zero

            vb = lx_a(lx_i + k + 1)
            vl = lx_a(lx_i + k + gxm)
            v  = lx_a(lx_i + k + gxm + 1)
            
            dvdx = (v - vl) / hx
            dvdy = (v - vb) / hy
            sqGrad = dvdx*dvdx + dvdy*dvdy
            avgWq = (2.0 * wq(i+1) + wq(i)) / 6.0d0
            avgWV = (wl(i+1)*vb + wl(i)*vl + wl(i+1)*v) / 3.0d0

            fu = avgWq * sqGrad - avgWV
            dvdx = dvdx * hy * avgWq
            dvdy = dvdy * hx * avgWq

            val1 = -dvdy - wl(i+1) * aread3
            call VecSetValuesLocal(G,1,k+1,val1,ADD_VALUES,info)
            
            val2 = -dvdx - wl(i) * aread3
            call VecSetValuesLocal(G,1,k+gxm,val2,ADD_VALUES,info)

            val3 = dvdx + dvdy - wl(i+1) * aread3
            call VecSetValuesLocal(G,1,k+gxm+1,val3,ADD_VALUES,info)

            fusum = fusum + fu
         enddo
      enddo

      floc = (fusum + flsum) * area
      call MPI_Allreduce(floc, f, 1, MPI_DOUBLE_PRECISION, MPI_SUM,               &
     &     MPI_COMM_WORLD, info)

      
      call VecRestoreArray(localX, lx_a, lx_i, info)
      call DARestoreLocalVector(da, localX, info)

      call VecAssemblyBegin(G, info)
      call VecAssemblyEnd(G, info)

      call PetscLogFlops((xe-xs) * (ye-ys) * 67 + 1,info)


      info = 0
      end


      subroutine WholeJBearHessian(daapp, da, X, H, dummy, info)
      implicit none
#include "jbearing3f.h"

      TAO_APPLICATION daapp
      DA              da
      Vec             X
      Mat             H
      PetscInt        dummy
      integer         info

      PetscInt        xs, xm, gxs, gxm
      PetscInt        xe, ys, ym
      PetscInt        gys, gym, ye
      PetscScalar     smallH(0:3, 0:3)

      PetscScalar     wql, wqu, dydx, dxdy
      PetscScalar     wqldydx, wqldxdy
      PetscScalar     wqudydx, wqudxdy
      PetscInt        i, j, ind(0:3)

      PetscScalar     wq_i, wq_ip1
      
      dydx = hy / hx
      dxdy = hx / hy

      call MatZeroEntries(H, info)

      call DAGetCorners(da, xs, ys, TAO_NULL_INTEGER, xm, ym,           &
     &     TAO_NULL_INTEGER, info)
      call DAGetGhostCorners(da, gxs, gys, TAO_NULL_INTEGER, gxm, gym,  &
     &     TAO_NULL_INTEGER, info)

      xe = gxs + gxm - 1
      ye = gys + gym - 1

      do j=ys,ye-1
         do i = xs, xe-1

            wq_i = wq(i)
            wq_ip1 = wq(i+1)

            wql = (2.0d0 * wq(i) + wq(i+1)) / 6.0d0
            wqu = (wq(i) + 2.0d0 * wq(i+1)) / 6.0d0
            

            wqldydx = wql * dydx
            wqldxdy = wql * dxdy
            wqudydx = wqu * dydx
            wqudxdy = wqu * dxdy

! Hessian contribution at 0,0
            smallH(0,0) = wqldxdy + wqldydx
            smallH(0,1) = -wqldydx
            smallH(1,0) = -wqldydx
            smallH(0,2) = -wqldxdy
            smallH(2,0) = -wqldxdy
            smallH(3,0) = 0.0d0
            smallH(0,3) = 0.0d0

! Hessian contribution at 1,0
            smallH(1,1) = wqldydx + wqudxdy
            smallH(1,2) = 0.0d0
            smallH(2,1) = 0.0d0
            smallH(1,3) = -wqudxdy
            smallH(3,1) = -wqudxdy

! Hessian contribution at 1,1
            smallH(2,2) = wqldxdy + wqudydx
            smallH(2,3) = -wqudydx
            smallH(3,2) = -wqudydx

! Hessian contribution at 1,1
            smallH(3,3) = wqudydx + wqudxdy

            ind(0) = (j-gys) * gxm + (i-gxs)
            ind(1) = ind(0) + 1
            ind(2) = ind(0) + gxm
            ind(3) = ind(2) + 1
            
            call MatSetValuesLocal(H,4,ind,4,ind,smallH,ADD_VALUES,info)
         enddo
      enddo

      call MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY, info)
      call MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY, info)
      
      call MatSetOption(H, MAT_SYMMETRIC, PETSC_TRUE,info)

      call PetscLogFlops((xe - xs) * (ye - ys) * 14 + 2,info)
      info = 0
      end

      
      

      
      
