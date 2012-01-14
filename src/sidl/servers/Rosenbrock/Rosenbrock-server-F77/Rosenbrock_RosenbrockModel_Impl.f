C       
C       File:          Rosenbrock_RosenbrockModel_Impl.f
C       Symbol:        Rosenbrock.RosenbrockModel-v0.0.2
C       Symbol Type:   class
C       Babel Version: 0.8.8
C       Description:   Server-side implementation for Rosenbrock.RosenbrockModel
C       
C       WARNING: Automatically generated; only changes within splicers preserved
C       
C       babel-version = 0.8.8
C       


C       
C       Symbol "Rosenbrock.RosenbrockModel" (version 0.0.2)
C       


C       DO-NOT-DELETE splicer.begin(_miscellaneous_code_start)
C       Insert extra code here...
C       DO-NOT-DELETE splicer.end(_miscellaneous_code_start)




C       
C       Class constructor called when the class is created.
C       

        subroutine Rosenbrock_RosenbrockModel__ctor_fi(self)
        implicit none
        integer*8 self

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel._ctor)
C       Insert the implementation here...
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel._ctor)
        end


C       
C       Class destructor called when the class is deleted.
C       

        subroutine Rosenbrock_RosenbrockModel__dtor_fi(self)
        implicit none
        integer*8 self

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel._dtor)
C       Insert the implementation here...
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel._dtor)
        end


C       
C       Method:  setNumberVariables[]
C       

        subroutine 
     &     Rosenbrock_RosenbrockModel_setNumberVariables_fi(self, n)
        implicit none
        integer*8 self
        integer*4 n

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.setNumberVariables)
        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/
        this_n = n;
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.setNumberVariables)
        end


C       
C       Method:  setAlpha[]
C       

        subroutine Rosenbrock_RosenbrockModel_setAlpha_fi(self, alpha)
        implicit none
        integer*8 self
        double precision alpha

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.setAlpha)
        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/
        this_alpha = alpha
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.setAlpha)
        end


C       
C       Optional --
C       initialize() will be called from the Optimization Solver 
C       before solving in order to set up any data.
C       

        subroutine Rosenbrock_RosenbrockModel_initialize_fi(self,
     &     retval)
        implicit none
        integer*8 self
        integer*4 retval

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.initialize)
        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/
        print *,'initializing rosenbrock (fortran server)'
        this_n = 2
        this_alpha = 99.0d0
        retval = 0
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.initialize)
        end


C       
C       Optional --
C       finalize() will be called from the Opimization Solver
C       when it through using the model (ie, when solver is deleted)
C       

        subroutine Rosenbrock_RosenbrockModel_finalize_fi(self, retval)
        implicit none
        integer*8 self
        integer*4 retval

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.finalize)
        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/

        retval = 0
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.finalize)
        end


C       
C       Required --
C       getNumberVariables()  lets the Optimization Solver know how much 
C       space to allocate for the various vectors and matrices.
C       

        subroutine 
     &     Rosenbrock_RosenbrockModel_getNumberVariables_fi(self,
     &     retval)
        implicit none
        integer*8 self
        integer*4 retval

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.getNumberVariables)
        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/

        retval = this_n
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.getNumberVariables)
        end


C       
C       Required --
C       evaluateObjectiveFunction(...) should be straightforward.
C       If you wish to only implement evaluateObjectiveAndGradient(), then
C       this function should allocate a dummy gradient vector and call
C       evaluateObjectiveAndGradient() explicitly.
C       

        subroutine 
     &     Rosenbrock_RosenbrockModel_evaluateObjectiveFunction_fi(self,
     &     x, f)
        implicit none
        integer*8 self
        integer*8 x
        double precision f

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.evaluateObjectiveFunction)
        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/

        double precision xraw(1), graw(1)
        common /arrays/ xraw, graw
        save /arrays/

        integer*4 xlower(1),xupper(1),xstride(1),xindex(1)
        integer*4 glower(1),gupper(1),gstride(1),gindex(1)


        call sidl_double__array_access_f(x, xraw, xlower, xupper, 
     &       xstride,xindex)


        if (xindex(1) .eq. 0 ) then
           print *,'Error in rosenbrock evaluation: array not aligned'
        endif

        call rawobjectivefunction(xraw(xindex(1)), f, this_n, 
     &       this_alpha)

C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.evaluateObjectiveFunction)
        end


C       
C       Required --
C       evaluateGradient(...) should be straightforward.
C       If you wish to only implement evaluateObjectiveAndGradient(), then
C       this function should call it explicitly. 
C       

        subroutine Rosenbrock_RosenbrockModel_evaluateGradient_fi(self,
     &     x, g)
        implicit none
        integer*8 self
        integer*8 x
        integer*8 g

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.evaluateGradient)

        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/
        
        double precision xraw(1), graw(1)
        common /arrays/ xraw, graw
        save /arrays/

        integer*4 xlower(1),xupper(1),xstride(1),xindex(1)
        integer*4 glower(1),gupper(1),gstride(1),gindex(1)


        call sidl_double__array_access_f(x,xraw,xlower,xupper,
     &       xstride,xindex)
        call sidl_double__array_access_f(g,graw,glower,gupper,
     &       gstride,gindex)



        if ((xindex(1) .eq. 0) .or. (gindex(1) .eq. 0)) then
           print *,'Error in rosenbrock evaluation: array not aligned'
        endif


        call rawgradient(xraw(xindex(1)), graw(gindex(1)), this_n,
     &       this_alpha)


C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.evaluateGradient)
        end


C       
C       Required --
C       evaluateObjectiveAndGradient(...) should be straightforward.
C       If you wish to implement the function and gradient evaluation routines
C       separately, then this method should do so explicitly.
C       

        subroutine 
     &     Rosenbrock_RosenbrockModel_evaluateObjectiveAndGradient_fi(
     &     self, x, f, g)
        implicit none
        integer*8 self
        integer*8 x
        double precision f
        integer*8 g

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.evaluateObjectiveAndGradient)

        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/
        
        double precision xraw(1), graw(1)
        common /arrays/ xraw, graw
        save /arrays/


        integer*4 xlower(1),xupper(1),xstride(1),xindex(1)
        integer*4 glower(1),gupper(1),gstride(1),gindex(1)


        
        call sidl_double__array_access_f(x,xraw,xlower,xupper,
     &       xstride,xindex)


        call sidl_double__array_access_f(g,graw,glower,gupper,
     &       gstride,gindex)

        if ((xindex(1) .eq. 0) .or. (gindex(1) .eq. 0)) then
           print *,'Error in rosenbrock evaluation: array not aligned'
        endif

        call rawobjectiveandgradient(xraw(xindex(1)),f,
     &          graw(gindex(1)), this_n,this_alpha)


C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.evaluateObjectiveAndGradient)
        end


C       
C       Optional --
C       evaluateHessian() is only necessary if using second derivative methods for
C       solving the model. 
C       

        subroutine Rosenbrock_RosenbrockModel_evaluateHessian_fi(self,
     &     x, H)
        implicit none
        integer*8 self
        integer*8 x
        integer*8 H

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.evaluateHessian)

        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/

        double precision xraw(1), graw(1)
        common /arrays/ xraw, graw

        integer*4 xlower(1), xupper(1), xstride(1), xindex(1)
        integer i

        call sidl_double__array_access_f(x, xraw, xlower, xupper,             &
     &     xstride, xindex)
        
        if (xindex(1) .eq. 0)  then
           print *,'Error in hessian evaluation: array not aligned'
        endif


        do i=0,this_n/2 - 1
           call sidl_double__array_set2_f(H, 2*i+1, 2*i+1, 
     &                                         2.0d0 * this_alpha)

           call sidl_double__array_set2_f(H, 2*i, 2*i,                         &
     &      -4.0d0 * this_alpha *                                              &
     &        ( xraw((2*i+1)*xstride(1) + xindex(1))                           &
     &                 -3.0d0 * xraw((2*i)*xstride(1) + xindex(1))**2)         &
     &      + 2.0d0)

           call sidl_double__array_set2_f(H,2*i+1, 2*i,                        &
     &        -4.0d0 * this_alpha * xraw((2*i)*xstride(1) + xindex(1)) )

           call sidl_double__array_set2_f(H,2*i, 2*i+1,                        &
     &        -4.0d0 * this_alpha * xraw((2*i)*xstride(1) + xindex(1)) )
                       
        enddo

C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.evaluateHessian)
        end


C       
C       Optional --
C       getVariableBounds() is only necessary if any of the variable bounds are 
C       not at (-inf, inf).    If a solution method
C       is selected that does not use variable bounds, then this function
C       will not be called.
C       

        subroutine Rosenbrock_RosenbrockModel_getVariableBounds_fi(self,
     &     xl, xu)
        implicit none
        integer*8 self
        integer*8 xl
        integer*8 xu

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.getVariableBounds)
        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/


C       Insert the implementation here...
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.getVariableBounds)
        end


C       
C       Optional --
C       initializeVariables() is called from the Optimization Solver to
C       sets the solution vector to an initial value.
C       

        subroutine 
     &     Rosenbrock_RosenbrockModel_initializeVariables_fi(self, x)
        implicit none
        integer*8 self
        integer*8 x

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.initializeVariables)

        integer this_n
        double precision this_alpha
        common /params/ this_alpha, this_n
        save /params/
        
        integer*4 i, lower,upper
        double precision val
        

        do i=0, this_n-1
           call Sidl_double__array_set1_f(x,i,0)
        enddo
C       Insert the implementation here...
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.initializeVariables)
        end


C       
C       Optional --
C       monitor() is called from TAO at every iteration.
C       If the application programmer wants to perform some kind of
C       visualization or other output throughout the solve, then this
C       method should be implemented.
C       

        subroutine Rosenbrock_RosenbrockModel_monitor_fi(self)
        implicit none
        integer*8 self

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.monitor)
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.monitor)
        end


C       
C       Optional --
C       checkConvergence() is called from TAO at every iteration.
C       If the application wishes to perform a convergence test, then
C       implement this method to set the flag to nonconverged, 0 for 
C       not converged.
C       The flag is an inout variable to avoid the need for implementing the
C       method if its utilization is not required (ie, 0 is passed in, so if
C       nothing is implemented, then 0 is passed back out).
C       

        subroutine Rosenbrock_RosenbrockModel_checkConvergence_fi(self,
     &     flag)
        implicit none
        integer*8 self
        integer*4 flag

C       DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.checkConvergence)
C       Insert the implementation here...
C       DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.checkConvergence)
        end


C       DO-NOT-DELETE splicer.begin(_miscellaneous_code_end)
      subroutine rawgradient(x,g,n,alpha)
      
      integer n
      double precision x(0:n-1), g(0:n-1), alpha

      double precision t1,t2
      integer i
      do i=0,n/2-1
         t1 = x(2*i + 1) - x(2*i) *x(2*i)
         t2 = 1.0d0 - x(2*i)
         g(2*i) = -4.0d0 * alpha * t1 * x(2*i) - 2.0d0*t2
         g(2*i+1) = 2.0d0 * alpha * t1
      enddo
      end


      subroutine rawobjectivefunction(x,f,n,alpha)

      integer n
      double precision x(0:n-1), f, alpha

      double precision t1,t2
      integer i
      f = 0
      do i=0,n/2-1
         t1 = x(2*i + 1) - x(2*i) *x(2*i)
         t2 = 1.0d0 - x(2*i)
         f = f + alpha*t1*t1 + t2*t2
      enddo
      end

      
      subroutine rawobjectiveandgradient(x,f,g,n,alpha)

      integer n
      double precision x(0:n-1),f,g(0:n-1),alpha
      

      double precision t1,t2
      integer i
      f = 0
      do i=0,n/2-1
         t1 = x(2*i+1) - x(2*i) * x(2*i)
         t2 = 1.0d0 - x(2*i)
         f = f + alpha*t1*t1 + t2*t2
         g(2*i) = -4.0d0 * alpha * t1 * x(2*i) - 2.0d0*t2
         g(2*i+1) = 2.0d0 * alpha * t1
      enddo
      end
         
C       DO-NOT-DELETE splicer.end(_miscellaneous_code_end)
