/*
 * File:          Rosenbrock_RosenbrockModel_Impl.c
 * Symbol:        Rosenbrock.RosenbrockModel-v0.0.2
 * Symbol Type:   class
 * Babel Version: 0.8.8
 * Description:   Server-side implementation for Rosenbrock.RosenbrockModel
 * 
 * WARNING: Automatically generated; only changes within splicers preserved
 * 
 * babel-version = 0.8.8
 */

/*
 * DEVELOPERS ARE EXPECTED TO PROVIDE IMPLEMENTATIONS
 * FOR THE FOLLOWING METHODS BETWEEN SPLICER PAIRS.
 */

/*
 * Symbol "Rosenbrock.RosenbrockModel" (version 0.0.2)
 */

#include "Rosenbrock_RosenbrockModel_Impl.h"

/* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel._includes) */
#include <stdio.h>
/* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel._includes) */

/*
 * Class constructor called when the class is created.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel__ctor"

void
impl_Rosenbrock_RosenbrockModel__ctor(
  Rosenbrock_RosenbrockModel self)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel._ctor) */

  /* Create private data structure */
  struct Rosenbrock_RosenbrockModel__data *dptr =
    (struct Rosenbrock_RosenbrockModel__data*)
        malloc(sizeof(struct Rosenbrock_RosenbrockModel__data));
  if (dptr) {
    dptr->n = 0;
    dptr->alpha = 0.0;
  }
  Rosenbrock_RosenbrockModel__set_data(self, dptr);

  return;
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel._ctor) */
}

/*
 * Class destructor called when the class is deleted.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel__dtor"

void
impl_Rosenbrock_RosenbrockModel__dtor(
  Rosenbrock_RosenbrockModel self)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel._dtor) */
  
  /* Free private data structure */
  struct Rosenbrock_RosenbrockModel__data *dptr =
            Rosenbrock_RosenbrockModel__get_data(self);

  if (dptr)
    free(dptr);
  else 
    fprintf(stderr,"Error: deleting nonexistant Rosenbrock_RosenbrockModel object\n");

  return;
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel._dtor) */
}

/*
 * Method:  setNumberVariables[]
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_setNumberVariables"

void
impl_Rosenbrock_RosenbrockModel_setNumberVariables(
  Rosenbrock_RosenbrockModel self, int32_t n)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.setNumberVariables) 
    */
  struct Rosenbrock_RosenbrockModel__data *dptr = 
    Rosenbrock_RosenbrockModel__get_data(self);

  dptr->n = n;

  Rosenbrock_RosenbrockModel__set_data(self,dptr);

  return;

  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.setNumberVariables) */
}

/*
 * Method:  setAlpha[]
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_setAlpha"

void
impl_Rosenbrock_RosenbrockModel_setAlpha(
  Rosenbrock_RosenbrockModel self, double alpha)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.setAlpha) */
  struct Rosenbrock_RosenbrockModel__data *dptr = 
    Rosenbrock_RosenbrockModel__get_data(self);

  dptr->alpha = alpha;

  Rosenbrock_RosenbrockModel__set_data(self,dptr);

  return;

  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.setAlpha) */
}

/*
 * Optional --
 * initialize() will be called from the Optimization Solver 
 * before solving in order to set up any data.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_initialize"

int32_t
impl_Rosenbrock_RosenbrockModel_initialize(
  Rosenbrock_RosenbrockModel self)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.initialize) */
  /* Get pointer to the private data structure */
  struct Rosenbrock_RosenbrockModel__data *dptr =
            Rosenbrock_RosenbrockModel__get_data(self);

  /* Set default data */
  dptr->n = 2;
  dptr->alpha = 99.0;

  Rosenbrock_RosenbrockModel__set_data(self,dptr);

  return 0;
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.initialize) */
}

/*
 * Optional --
 * finalize() will be called from the Opimization Solver
 * when it through using the model (ie, when solver is deleted)
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_finalize"

int32_t
impl_Rosenbrock_RosenbrockModel_finalize(
  Rosenbrock_RosenbrockModel self)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.finalize) */
  return 0;
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.finalize) */
}

/*
 * Required --
 * getNumberVariables()  lets the Optimization Solver know how much 
 * space to allocate for the various vectors and matrices.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_getNumberVariables"

int32_t
impl_Rosenbrock_RosenbrockModel_getNumberVariables(
  Rosenbrock_RosenbrockModel self)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.getNumberVariables) 
    */
  struct Rosenbrock_RosenbrockModel__data *dptr =
            Rosenbrock_RosenbrockModel__get_data(self);
  return dptr->n;
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.getNumberVariables) */
}

/*
 * Required --
 * evaluateObjectiveFunction(...) should be straightforward.
 * If you wish to only implement evaluateObjectiveAndGradient(), then
 * this function should allocate a dummy gradient vector and call
 * evaluateObjectiveAndGradient() explicitly.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_evaluateObjectiveFunction"

void
impl_Rosenbrock_RosenbrockModel_evaluateObjectiveFunction(
  Rosenbrock_RosenbrockModel self, struct SIDL_double__array* x, double* f)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.evaluateObjectiveFunction) */
  struct Rosenbrock_RosenbrockModel__data *dptr =
            Rosenbrock_RosenbrockModel__get_data(self);

  double t1,t2;
  int i;
  double *rawx = x->d_firstElement;
  
  *f = 0;
  for (i=0;i<dptr->n/2; i++) {
    t1 = rawx[2*i+1] - rawx[2*i]*rawx[2*i];
    t2 = 1-rawx[2*i];
    *f += dptr->alpha * t1*t1 + t2*t2;
  }
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.evaluateObjectiveFunction) */
}

/*
 * Required --
 * evaluateGradient(...) should be straightforward.
 * If you wish to only implement evaluateObjectiveAndGradient(), then
 * this function should call it explicitly. 
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_evaluateGradient"

void
impl_Rosenbrock_RosenbrockModel_evaluateGradient(
  Rosenbrock_RosenbrockModel self, struct SIDL_double__array* x,
    struct SIDL_double__array* g)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.evaluateGradient) */
  struct Rosenbrock_RosenbrockModel__data *dptr =
            Rosenbrock_RosenbrockModel__get_data(self);

  double t1,t2;
  int i;
  double *rawx = x->d_firstElement;
  double *rawg = g->d_firstElement;
  
  for (i=0;i<dptr->n/2; i++) {
    t1 = rawx[2*i+1] - rawx[2*i]*rawx[2*i];
    t2 = 1-rawx[2*i];
    rawg[2*i] = -4.0*dptr->alpha * t1* rawx[2*i] - 2.0*t2;
    rawg[2*i+1] = 2.0*dptr->alpha * t1;
  }
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.evaluateGradient) */
}

/*
 * Required --
 * evaluateObjectiveAndGradient(...) should be straightforward.
 * If you wish to implement the function and gradient evaluation routines
 * separately, then this method should do so explicitly.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_evaluateObjectiveAndGradient"

void
impl_Rosenbrock_RosenbrockModel_evaluateObjectiveAndGradient(
  Rosenbrock_RosenbrockModel self, struct SIDL_double__array* x, double* f,
    struct SIDL_double__array* g)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.evaluateObjectiveAndGradient) */
  struct Rosenbrock_RosenbrockModel__data *dptr =
    Rosenbrock_RosenbrockModel__get_data(self);

  double t1, t2;
  int i;
  double *rawx = x->d_firstElement;
  double *rawg = g->d_firstElement;

  *f = 0;

  for (i=0; i<dptr->n/2; i++) {
    t1 = rawx[2*i+1] - rawx[2*i]*rawx[2*i];
    t2 = 1-rawx[2*i];
    *f += dptr->alpha * t1*t1 + t2*t2;
    rawg[2*i] = -4.0*dptr->alpha * t1* rawx[2*i] - 2.0*t2;
    rawg[2*i+1] = 2.0*dptr->alpha * t1;
  }

  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.evaluateObjectiveAndGradient) */
}

/*
 * Optional --
 * evaluateHessian() is only necessary if using second derivative methods for
 * solving the model. 
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_evaluateHessian"

void
impl_Rosenbrock_RosenbrockModel_evaluateHessian(
  Rosenbrock_RosenbrockModel self, struct SIDL_double__array* x,
    struct SIDL_double__array* H)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.evaluateHessian) */
  struct Rosenbrock_RosenbrockModel__data *dptr =
            Rosenbrock_RosenbrockModel__get_data(self);
  int i;
  double xeven, xodd;
  
  for (i=0;i<dptr->n/2;i++) {
    xeven = sidlArrayElem1(x,2*i);
    xodd = sidlArrayElem1(x,2*i+1);
    
    sidlArrayElem2(H,2*i+1, 2*i+1) = 2*dptr->alpha;
    sidlArrayElem2(H,2*i, 2*i) = 2.0 - 4.0 * dptr->alpha * (xodd - 3.0*xeven*xeven);
    sidlArrayElem2(H,2*i+1,2*i) = sidlArrayElem2(H,2*i,2*i+1) = 
      -4.0 * dptr->alpha * xeven;
  }

  return;
    
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.evaluateHessian) */
}

/*
 * Optional --
 * getVariableBounds() is only necessary if any of the variable bounds are 
 * not at (-inf, inf).    If a solution method
 * is selected that does not use variable bounds, then this function
 * will not be called.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_getVariableBounds"

void
impl_Rosenbrock_RosenbrockModel_getVariableBounds(
  Rosenbrock_RosenbrockModel self, struct SIDL_double__array* xl,
    struct SIDL_double__array* xu)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.getVariableBounds) 
    */

  return;
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.getVariableBounds) */
}

/*
 * Optional --
 * initializeVariables() is called from the Optimization Solver to
 * sets the solution vector to an initial value.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_initializeVariables"

void
impl_Rosenbrock_RosenbrockModel_initializeVariables(
  Rosenbrock_RosenbrockModel self, struct SIDL_double__array* x)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.initializeVariables) */
  int i;
  struct Rosenbrock_RosenbrockModel__data *dptr =
            Rosenbrock_RosenbrockModel__get_data(self);
  
  for (i=0;i<dptr->n;i++)
    sidlArrayElem1(x,i) = 0.0;
  return;
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.initializeVariables) */
}

/*
 * Optional --
 * monitor() is called from TAO at every iteration.
 * If the application programmer wants to perform some kind of
 * visualization or other output throughout the solve, then this
 * method should be implemented.
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_monitor"

void
impl_Rosenbrock_RosenbrockModel_monitor(
  Rosenbrock_RosenbrockModel self)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.monitor) */
  return;
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.monitor) */
}

/*
 * Optional --
 * checkConvergence() is called from TAO at every iteration.
 * If the application wishes to perform a convergence test, then
 * implement this method to set the flag to nonconverged, 0 for 
 * not converged.
 * The flag is an inout variable to avoid the need for implementing the
 * method if its utilization is not required (ie, 0 is passed in, so if
 * nothing is implemented, then 0 is passed back out).
 */

#undef __FUNC__
#define __FUNC__ "impl_Rosenbrock_RosenbrockModel_checkConvergence"

void
impl_Rosenbrock_RosenbrockModel_checkConvergence(
  Rosenbrock_RosenbrockModel self, int32_t* flag)
{
  /* DO-NOT-DELETE splicer.begin(Rosenbrock.RosenbrockModel.checkConvergence) */
  /* Insert the implementation of the checkConvergence method here... */
  /* DO-NOT-DELETE splicer.end(Rosenbrock.RosenbrockModel.checkConvergence) */
}
