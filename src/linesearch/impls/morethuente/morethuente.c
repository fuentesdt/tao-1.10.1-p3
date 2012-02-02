/*$Id$*/

#include "morethuente.h"  /*I "tao_solver.h" I*/
#include <stdlib.h>

/*
     The subroutine mcstep is taken from the work of Jorge Nocedal.
     this is a variant of More' and Thuente's routine.

     subroutine mcstep

     the purpose of mcstep is to compute a safeguarded step for
     a linesearch and to update an interval of uncertainty for
     a minimizer of the function.

     the parameter stx contains the step with the least function
     value. the parameter stp contains the current step. it is
     assumed that the derivative at stx is negative in the
     direction of the step. if bracket is set true then a
     minimizer has been bracketed in an interval of uncertainty
     with endpoints stx and sty.

     the subroutine statement is

     subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,bracket,
                       stpmin,stpmax,info)

     where

       stx, fx, and dx are variables which specify the step,
         the function, and the derivative at the best step obtained
         so far. The derivative must be negative in the direction
         of the step, that is, dx and stp-stx must have opposite
         signs. On output these parameters are updated appropriately.

       sty, fy, and dy are variables which specify the step,
         the function, and the derivative at the other endpoint of
         the interval of uncertainty. On output these parameters are
         updated appropriately.

       stp, fp, and dp are variables which specify the step,
         the function, and the derivative at the current step.
         If bracket is set true then on input stp must be
         between stx and sty. On output stp is set to the new step.

       bracket is a logical variable which specifies if a minimizer
         has been bracketed.  If the minimizer has not been bracketed
         then on input bracket must be set false.  If the minimizer
         is bracketed then on output bracket is set true.

       stpmin and stpmax are input variables which specify lower
         and upper bounds for the step.

       info is an integer output variable set as follows:
         if info = 1,2,3,4,5, then the step has been computed
         according to one of the five cases below. otherwise
         info = 0, and this indicates improper input parameters.

     subprograms called

       fortran-supplied ... abs,max,min,sqrt

     argonne national laboratory. minpack project. june 1983
     jorge j. more', david j. thuente

*/

#undef __FUNCT__  
#define __FUNCT__ "TaoStep_LineSearch"
static int TaoStep_LineSearch(TAO_SOLVER tao,
                              double *stx, double *fx, double *dx,
		              double *sty, double *fy, double *dy,
			      double *stp, double *fp, double *dp)
{
  TAO_LINESEARCH *neP = (TAO_LINESEARCH *) tao->linectx;
  double gamma1, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  TaoInt bound;

  TaoFunctionBegin;

  // Check the input parameters for errors
  neP->infoc = 0;
  if (neP->bracket && (*stp <= TaoMin(*stx,*sty) || (*stp >= TaoMax(*stx,*sty)))) SETERRQ(1,"bad stp in bracket");
  if (*dx * (*stp-*stx) >= 0.0) SETERRQ(1,"dx * (stp-stx) >= 0.0");
  if (neP->stepmax < neP->stepmin) SETERRQ(1,"stepmax > stepmin");

  // Determine if the derivatives have opposite sign */
  sgnd = *dp * (*dx / TaoAbsDouble(*dx));

  if (*fp > *fx) {
    // Case 1: a higher function value.
    // The minimum is bracketed. If the cubic step is closer
    // to stx than the quadratic step, the cubic step is taken,
    // else the average of the cubic and quadratic steps is taken.

    neP->infoc = 1;
    bound = 1;
    theta = 3 * (*fx - *fp) / (*stp - *stx) + *dx + *dp;
    s = TaoMax(TaoAbsDouble(theta),TaoAbsDouble(*dx));
    s = TaoMax(s,TaoAbsDouble(*dp));
    gamma1 = s*sqrt(pow(theta/s,2.0) - (*dx/s)*(*dp/s));
    if (*stp < *stx) gamma1 = -gamma1;
    /* Can p be 0?  Check */
    p = (gamma1 - *dx) + theta;
    q = ((gamma1 - *dx) + gamma1) + *dp;
    r = p/q;
    stpc = *stx + r*(*stp - *stx);
    stpq = *stx + ((*dx/((*fx-*fp)/(*stp-*stx)+*dx))*0.5) * (*stp - *stx);

    if (TaoAbsDouble(stpc-*stx) < TaoAbsDouble(stpq-*stx)) {
      stpf = stpc;
    } 
    else {
      stpf = stpc + 0.5*(stpq - stpc);
    }
    neP->bracket = 1;
  }
  else if (sgnd < 0.0) {
    // Case 2: A lower function value and derivatives of
    // opposite sign. The minimum is bracketed. If the cubic
    // step is closer to stx than the quadratic (secant) step,
    // the cubic step is taken, else the quadratic step is taken.

    neP->infoc = 2;
    bound = 0;
    theta = 3*(*fx - *fp)/(*stp - *stx) + *dx + *dp;
    s = TaoMax(TaoAbsDouble(theta),TaoAbsDouble(*dx));
    s = TaoMax(s,TaoAbsDouble(*dp));
    gamma1 = s*sqrt(pow(theta/s,2.0) - (*dx/s)*(*dp/s));
    if (*stp > *stx) gamma1 = -gamma1;
    p = (gamma1 - *dp) + theta;
    q = ((gamma1 - *dp) + gamma1) + *dx;
    r = p/q;
    stpc = *stp + r*(*stx - *stp);
    stpq = *stp + (*dp/(*dp-*dx))*(*stx - *stp);

    if (TaoAbsDouble(stpc-*stp) > TaoAbsDouble(stpq-*stp)) {
      stpf = stpc;
    }
    else {
      stpf = stpq;
    }
    neP->bracket = 1;
  }
  else if (TaoAbsDouble(*dp) < TaoAbsDouble(*dx)) {
    // Case 3: A lower function value, derivatives of the
    // same sign, and the magnitude of the derivative decreases.
    // The cubic step is only used if the cubic tends to infinity
    // in the direction of the step or if the minimum of the cubic
    // is beyond stp. Otherwise the cubic step is defined to be
    // either stepmin or stepmax. The quadratic (secant) step is also
    // computed and if the minimum is bracketed then the the step
    // closest to stx is taken, else the step farthest away is taken.

    neP->infoc = 3;
    bound = 1;
    theta = 3*(*fx - *fp)/(*stp - *stx) + *dx + *dp;
    s = TaoMax(TaoAbsDouble(theta),TaoAbsDouble(*dx));
    s = TaoMax(s,TaoAbsDouble(*dp));

    // The case gamma1 = 0 only arises if the cubic does not tend
    // to infinity in the direction of the step.
    gamma1 = s*sqrt(TaoMax(0.0,pow(theta/s,2.0) - (*dx/s)*(*dp/s)));
    if (*stp > *stx) gamma1 = -gamma1;
    p = (gamma1 - *dp) + theta;
    q = (gamma1 + (*dx - *dp)) + gamma1;
    r = p/q;
    if (r < 0.0 && gamma1 != 0.0) stpc = *stp + r*(*stx - *stp);
    else if (*stp > *stx)        stpc = neP->stepmax;
    else                         stpc = neP->stepmin;
    stpq = *stp + (*dp/(*dp-*dx)) * (*stx - *stp);

    if (neP->bracket) {
      if (TaoAbsDouble(*stp-stpc) < TaoAbsDouble(*stp-stpq)) {
	stpf = stpc;
      } 
      else {
	stpf = stpq;
      }
    }
    else {
      if (TaoAbsDouble(*stp-stpc) > TaoAbsDouble(*stp-stpq)) {
	stpf = stpc;
      }
      else {
	stpf = stpq;
      }
    }
  }
  else {
    // Case 4: A lower function value, derivatives of the
    // same sign, and the magnitude of the derivative does
    // not decrease. If the minimum is not bracketed, the step
    // is either stpmin or stpmax, else the cubic step is taken.

    neP->infoc = 4;
    bound = 0;
    if (neP->bracket) {
      theta = 3*(*fp - *fy)/(*sty - *stp) + *dy + *dp;
      s = TaoMax(TaoAbsDouble(theta),TaoAbsDouble(*dy));
      s = TaoMax(s,TaoAbsDouble(*dp));
      gamma1 = s*sqrt(pow(theta/s,2.0) - (*dy/s)*(*dp/s));
      if (*stp > *sty) gamma1 = -gamma1;
      p = (gamma1 - *dp) + theta;
      q = ((gamma1 - *dp) + gamma1) + *dy;
      r = p/q;
      stpc = *stp + r*(*sty - *stp);
      stpq = *stp + (*dp/(*dp-*dx)) * (*stx - *stp);

      stpf = stpc;
    } 
    else if (*stp > *stx) {
      stpf = neP->stepmax;
    } 
    else {
      stpf = neP->stepmin;
    }
  }
  
  // Update the interval of uncertainty.  This update does not
  // depend on the new step or the case analysis above.

  if (*fp > *fx) {
    *sty = *stp;
    *fy = *fp;
    *dy = *dp;
  } 
  else {
    if (sgnd < 0.0) {
      *sty = *stx;
      *fy = *fx;
      *dy = *dx;
    }
    *stx = *stp;
    *fx = *fp;
    *dx = *dp;
  }
  
  // Compute the new step and safeguard it.
  stpf = TaoMin(neP->stepmax,stpf);
  stpf = TaoMax(neP->stepmin,stpf);
  *stp = stpf;
  if (neP->bracket && bound) {
    if (*sty > *stx) {
      *stp = TaoMin(*stx+0.66*(*sty-*stx),*stp);
    }
    else {
      *stp = TaoMax(*stx+0.66*(*sty-*stx),*stp);
    }
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_LineSearch"
/* @ 
   TaoApply_LineSearch - This routine performs a line search algorithm.

   Input Parameters:
+  tao - TAO_SOLVER context
.  X - current iterate (on output X contains new iterate, X + step*S)
.  S - search direction
.  f - objective function evaluated at X
.  G - gradient evaluated at X
.  W - work vector
-  step - initial estimate of step length

   Output parameters:
+  f - objective function evaluated at new iterate, X + step*S
.  G - gradient evaluated at new iterate, X + step*S
.  X - new iterate
.  gnorm - 2-norm of G
-  step - final step length


   Info is set to one of:
+   0 - the line search succeeds; the sufficient decrease
   condition and the directional derivative condition hold

   negative number if an input parameter is invalid
.   -1 -  step < 0 
.   -2 -  ftol < 0 
.   -3 -  rtol < 0 
.   -4 -  gtol < 0 
.   -5 -  stepmin < 0 
.   -6 -  stepmax < stepmin
-   -7 -  maxfev < 0

   positive number > 1 if the line search otherwise terminates
+    2 -  Relative width of the interval of uncertainty is 
         at most rtol.
.    3 -  Maximum number of function evaluations (maxfev) has 
         been reached.
.    4 -  Step is at the lower bound, stepmin.
.    5 -  Step is at the upper bound, stepmax.
.    6 -  Rounding errors may prevent further progress. 
         There may not be a step that satisfies the 
         sufficient decrease and curvature conditions.  
         Tolerances may be too small.
-    7 -  Search direction is not a descent direction.

   Notes:
   This algorithm is taken from More' and Thuente, "Line search algorithms
   with guaranteed sufficient decrease", Argonne National Laboratory, 
   Technical Report MCS-P330-1092.

   Notes:
   This routine is used within the following TAO unconstrained minimization 
   solvers: Newton linesearch (tao_nls), limited memory variable metric
   (tao_lmvm), and conjugate gradient (tao_cg).

   Level: advanced

.keywords: TAO_SOLVER, linesearch
@ */
static int TaoApply_LineSearch(TAO_SOLVER tao,TaoVec* X,TaoVec* G,TaoVec* S,TaoVec* W,double *f, double *f_full,
                        double *step,TaoInt *info2,void*ctx)
{
  TAO_LINESEARCH *neP = (TAO_LINESEARCH *) tao->linectx;
  double    xtrapf = 4.0;
  double    finit, width, width1, dginit, fm, fxm, fym, dgm, dgxm, dgym;
  double    dgx, dgy, dg, fx, fy, stx, sty, dgtest, ftest1=0.0;
  int       info;
  TaoInt  i, stage1;
  int    need_gradient=0;

#if defined(PETSC_USE_COMPLEX)
  PetscScalar    cdginit, cdg, cstep = 0.0;
#endif

  TaoFunctionBegin;
  /* neP->stepmin - lower bound for step */
  /* neP->stepmax - upper bound for step */
  /* neP->rtol 	  - relative tolerance for an acceptable step */
  /* neP->ftol 	  - tolerance for sufficient decrease condition */
  /* neP->gtol 	  - tolerance for curvature condition */
  /* neP->nfev 	  - number of function evaluations */
  /* neP->maxfev  - maximum number of function evaluations */

  *info2 = 0;

  /* Check input parameters for errors */
  if (*step < 0.0) {
    info = PetscInfo1(tao, "TaoApply_LineSearch: Line search error: step (%g) < 0\n",*step); CHKERRQ(info);
    *info2 = -1;
  } 
  
  if (neP->ftol < 0.0) {
    info = PetscInfo1(tao, "TaoApply_LineSearch: Line search error: ftol (%g) < 0\n",neP->ftol); CHKERRQ(info);
    *info2 = -2;
  } 
  
  if (neP->rtol < 0.0) {
    info = PetscInfo1(tao, "TaoApply_LineSearch: Line search error: rtol (%g) < 0\n",neP->rtol); CHKERRQ(info);
    *info2 = -3;
  } 

  if (neP->gtol < 0.0) {
    info = PetscInfo1(tao, "TaoApply_LineSearch: Line search error: gtol (%g) < 0\n",neP->gtol); CHKERRQ(info);
    *info2 = -4;
  } 

  if (neP->stepmin < 0.0) {
    info = PetscInfo1(tao, "TaoApply_LineSearch: Line search error: stepmin (%g) < 0\n",neP->stepmin); CHKERRQ(info);
    *info2 = -5;
  } 

  if (neP->stepmax < neP->stepmin) {
    info = PetscInfo2(tao, "TaoApply_LineSearch: Line search error: stepmax (%g) < stepmin (%g)\n", neP->stepmax,neP->stepmin); CHKERRQ(info);
    *info2 = -6;
  }
  
  if (neP->maxfev < 0) {
    info = PetscInfo1(tao, "TaoApply_LineSearch: Line search error: maxfev (%d) < 0\n",neP->maxfev); CHKERRQ(info);
    *info2 = -7;
  }

  if (TaoInfOrNaN(*f)) {
    info = PetscInfo1(tao, "TaoApply_LineSearch: Line search error: function (%g) inf or nan\n", *f); CHKERRQ(info);
    *info2 = -8;
  }

  /* Check that search direction is a descent direction */

#if defined(PETSC_USE_COMPLEX)
  info = G->Dot(S,&cdginit);CHKERRQ(info); dginit = TaoReal(cdginit);
#else
  info = G->Dot(S,&dginit);CHKERRQ(info);
#endif

  if (TaoInfOrNaN(dginit)) {
    info = PetscInfo1(tao,"TaoApply_LineSearch: Line search error: dginit (%g) inf or nan\n", dginit); CHKERRQ(info);
    *info2 = -9;
  }

  if (dginit >= 0.0) {
    info = PetscInfo(tao,"TaoApply_LineSearch:Search direction not a descent direction\n"); CHKERRQ(info);
    *info2 = -10;
  }

  if (*info2) {
    TaoFunctionReturn(0);
  }

  /* Initialization */
  neP->bracket = 0;
  *info2 = 0;
  stage1 = 1;
  finit = *f;
  dgtest = neP->ftol * dginit;
  width = neP->stepmax - neP->stepmin;
  width1 = width * 2.0;
  info = W->CopyFrom(X);CHKERRQ(info);
  /* Variable dictionary:  
     stx, fx, dgx - the step, function, and derivative at the best step
     sty, fy, dgy - the step, function, and derivative at the other endpoint 
                   of the interval of uncertainty
     step, f, dg - the step, function, and derivative at the current step */

  stx = 0.0;
  fx  = finit;
  dgx = dginit;
  sty = 0.0;
  fy  = finit;
  dgy = dginit;
 
  neP->nfev = 0;
  tao->new_search=TAO_TRUE;
  for (i=0; i< neP->maxfev; i++) {
    /* Set min and max steps to correspond to the interval of uncertainty */
    if (neP->bracket) {
      neP->stepmin = TaoMin(stx,sty); 
      neP->stepmax = TaoMax(stx,sty); 
    } 
    else {
      neP->stepmin = stx;
      neP->stepmax = *step + xtrapf * (*step - stx);
    }

    /* Force the step to be within the bounds */
    *step = TaoMax(*step,neP->stepmin);
    *step = TaoMin(*step,neP->stepmax);
    
    /* If an unusual termination is to occur, then let step be the lowest
       point obtained thus far */
    if (((neP->bracket) && (*step <= neP->stepmin || *step >= neP->stepmax)) ||
        ((neP->bracket) && (neP->stepmax - neP->stepmin <= neP->rtol * neP->stepmax)) ||
        (neP->nfev >= neP->maxfev - 1) || (neP->infoc == 0)) {
      *step = stx;
    }

#if defined(PETSC_USE_COMPLEX)
    cstep = *step;
    info = W->Waxpby(cstep,S,1.0,X);CHKERRQ(info);
#else
    info = W->Waxpby(*step,S,1.0,X);CHKERRQ(info); 	/* W = X + step*S */
#endif
    tao->current_step=*step;
    if (tao->MeritFunctionGTSApply) {
      info = TaoComputeMeritFunctionGTS(tao,W,S,f,&dg); CHKERRQ(info);
      need_gradient = 1;
    } else {
      info = TaoComputeMeritFunctionGradient(tao,W,f,G);CHKERRQ(info);
      tao->new_search=TAO_FALSE;
      neP->nfev++;
#if defined(PETSC_USE_COMPLEX)
      info = G->Dot(S,&cdg);CHKERRQ(info); dg = TaoReal(cdg);
#else
      info = G->Dot(S,&dg);CHKERRQ(info);	        /* dg = G^T S */
#endif
    }
    if (0 == i) {
      *f_full = *f;
    }
      
    if (TaoInfOrNaN(*f) || TaoInfOrNaN(dg)) {
      // User provided compute function generated Not-a-Number, assume 
      // domain violation and set function value and directional
      // derivative to infinity.
      *f = TAO_INFINITY;
      dg = TAO_INFINITY;
    }

    ftest1 = finit + *step * dgtest;

    /* Convergence testing */
    if (((*f - ftest1 <= 1.0e-10 * fabs(finit)) && 
         (TaoAbsDouble(dg) + neP->gtol*dginit <= 0.0))) {
      info = PetscInfo(tao, "TaoApply_LineSearch:Line search success: Sufficient decrease and directional deriv conditions hold\n"); CHKERRQ(info);
      *info2 = 0; 
      break;
    }

    /* Checks for bad cases */
    if (((neP->bracket) && (*step <= neP->stepmin||*step >= neP->stepmax)) || (!neP->infoc)) {
      info = PetscInfo(tao,"TaoApply_LineSearch:Rounding errors may prevent further progress.  May not be a step satisfying\n"); CHKERRQ(info);
      info = PetscInfo(tao,"TaoApply_LineSearch:sufficient decrease and curvature conditions. Tolerances may be too small.\n"); CHKERRQ(info);
      *info2 = 6; break;
    }
    if ((*step == neP->stepmax) && (*f <= ftest1) && (dg <= dgtest)) {
      info = PetscInfo1(tao,"TaoApply_LineSearch:Step is at the upper bound, stepmax (%g)\n",neP->stepmax); CHKERRQ(info);
      *info2 = 5; break;
    }
    if ((*step == neP->stepmin) && (*f >= ftest1) && (dg >= dgtest)) {
      info = PetscInfo1(tao,"TaoApply_LineSearch:Step is at the lower bound, stepmin (%g)\n",neP->stepmin); CHKERRQ(info);
      *info2 = 4; break;
    }
    if (neP->nfev >= neP->maxfev) {
      info = PetscInfo2(tao,"TaoApply_LineSearch:Number of line search function evals (%d) > maximum (%d)\n",neP->nfev,neP->maxfev); CHKERRQ(info);
      *info2 = 3; break;
    }
    if ((neP->bracket) && (neP->stepmax - neP->stepmin <= neP->rtol*neP->stepmax)){
      info = PetscInfo1(tao,"TaoApply_LineSearch:Relative width of interval of uncertainty is at most rtol (%g)\n",neP->rtol); CHKERRQ(info);
      *info2 = 2; break;
    }

    /* In the first stage, we seek a step for which the modified function
       has a nonpositive value and nonnegative derivative */
    if ((stage1) && (*f <= ftest1) && (dg >= dginit * TaoMin(neP->ftol, neP->gtol))) {
      stage1 = 0;
    }

    /* A modified function is used to predict the step only if we
       have not obtained a step for which the modified function has a 
       nonpositive function value and nonnegative derivative, and if a
       lower function value has been obtained but the decrease is not
       sufficient */

    if ((stage1) && (*f <= fx) && (*f > ftest1)) {
      fm   = *f - *step * dgtest;	/* Define modified function */
      fxm  = fx - stx * dgtest;	        /* and derivatives */
      fym  = fy - sty * dgtest;
      dgm  = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;
      
      /* Update the interval of uncertainty and compute the new step */
      info = TaoStep_LineSearch(tao,&stx,&fxm,&dgxm,&sty,&fym,&dgym,step,&fm,&dgm);CHKERRQ(info);
      
      fx  = fxm + stx * dgtest;	/* Reset the function and */
      fy  = fym + sty * dgtest;	/* gradient values */
      dgx = dgxm + dgtest; 
      dgy = dgym + dgtest; 
    } 
    else {
      /* Update the interval of uncertainty and compute the new step */
      info = TaoStep_LineSearch(tao,&stx,&fx,&dgx,&sty,&fy,&dgy,step,f,&dg);CHKERRQ(info);
    }
    
    /* Force a sufficient decrease in the interval of uncertainty */
    if (neP->bracket) {
      if (TaoAbsDouble(sty - stx) >= 0.66 * width1) *step = stx + 0.5*(sty - stx);
      width1 = width;
      width = TaoAbsDouble(sty - stx);
    }
  }
  
  /* Finish computations */
  info = X->CopyFrom(W); CHKERRQ(info);
  if (need_gradient) {
    info = TaoComputeMeritGradient(tao,X,G); CHKERRQ(info);
  }
  info = PetscInfo2(tao,"TaoApply_LineSearch:%d function evals in line search, step = %10.4f\n",neP->nfev,*step); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_BoundLineSearch"
/* @ 
   TaoApply_BoundLineSearch - This routine performs a line search algorithm.

   Input Parameters:
+  tao - TAO_SOLVER context
.  X - current iterate (on output X contains new iterate, X + step*S)
.  S - search direction
.  f - objective function evaluated at X
.  G - gradient evaluated at X
.  W - work vector
-  step - initial estimate of step length

   Output parameters:
+  f - objective function evaluated at new iterate, X + step*S
.  G - gradient evaluated at new iterate, X + step*S
.  X - new iterate
.  gnorm - 2-norm of G
-  step - final step length


   Info is set to one of:
+   0 - the line search succeeds; the sufficient decrease
   condition and the directional derivative condition hold

   negative number if an input parameter is invalid
.   -1 -  step < 0 
.   -2 -  ftol < 0 
.   -3 -  rtol < 0 
.   -4 -  gtol < 0 
.   -5 -  stepmin < 0 
.   -6 -  stepmax < stepmin
-   -7 -  maxfev < 0

   positive number > 1 if the line search otherwise terminates
+    2 -  Relative width of the interval of uncertainty is 
         at most rtol.
.    3 -  Maximum number of function evaluations (maxfev) has 
         been reached.
.    4 -  Step is at the lower bound, stepmin.
.    5 -  Step is at the upper bound, stepmax.
.    6 -  Rounding errors may prevent further progress. 
         There may not be a step that satisfies the 
         sufficient decrease and curvature conditions.  
         Tolerances may be too small.
-    7 -  Search direction is not a descent direction.

   Notes:
   This algorithm is is a modification of the algorithm by More' and
   Thuente.  The modifications concern bounds.  This algorithm steps
   in the direction passed into this routine.  This point get
   projected back into the feasible set.  In the context of bound
   constrained optimization, there may not be a point in the piecewise
   linear path that satisfies the Wolfe conditions.  When the active
   set is changing, decrease in the objective function may be
   sufficient to terminate this line search.

   Note: Much of this coded is identical to the More' Thuente line search.
   Variations to the code are commented.

   Notes:
   This routine is used within the following TAO bound constrained minimization
   solvers: Newton linesearch (tao_tron) and limited memory variable metric
   (tao_blmvm).

   Level: advanced

.keywords: TAO_SOLVER, linesearch
@ */
static int TaoApply_BoundLineSearch(TAO_SOLVER tao,TaoVec* X,TaoVec* G,TaoVec* S,TaoVec* W,double *f, double *f_full,
                        double *step,TaoInt *info2,void*ctx)
{
  TAO_LINESEARCH *neP = (TAO_LINESEARCH *) tao->linectx;
  TaoVec *XL,*XU;
  double    xtrapf = 4.0;
  double    finit, width, width1, dginit, fm, fxm, fym, dgm, dgxm, dgym;
  double    dgx, dgy, dg, fx, fy, stx, sty, dgtest, ftest1=0.0, ftest2=0.0;
  double    bstepmin1, bstepmin2, bstepmax;
  double    dg1, dg2;
  int       info;
  int       need_gradient=0;
  TaoInt    i, stage1;

#if defined(PETSC_USE_COMPLEX)
  PetscScalar    cdginit, cdg, cstep = 0.0;
#endif

  TaoFunctionBegin;
  /* neP->stepmin - lower bound for step */
  /* neP->stepmax - upper bound for step */
  /* neP->rtol 	  - relative tolerance for an acceptable step */
  /* neP->ftol 	  - tolerance for sufficient decrease condition */
  /* neP->gtol 	  - tolerance for curvature condition */
  /* neP->nfev 	  - number of function evaluations */
  /* neP->maxfev  - maximum number of function evaluations */

  /* Check input parameters for errors */
  if (*step < 0.0) {
    info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Line search error: step (%g) < 0\n",*step); CHKERRQ(info);
    *info2 = -1; TaoFunctionReturn(0);
  } 
  else if (neP->ftol < 0.0) {
    info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Line search error: ftol (%g) < 0\n",neP->ftol); CHKERRQ(info);
    *info2 = -2; TaoFunctionReturn(0);
  } 
  else if (neP->rtol < 0.0) {
    info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Line search error: rtol (%g) < 0\n",neP->rtol); CHKERRQ(info);
    *info2 = -3; TaoFunctionReturn(0);
  } 
  else if (neP->gtol < 0.0) {
    info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Line search error: gtol (%g) < 0\n",neP->gtol); CHKERRQ(info);
    *info2 = -4; TaoFunctionReturn(0);
  } 
  else if (neP->stepmin < 0.0) {
    info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Line search error: stepmin (%g) < 0\n",neP->stepmin); CHKERRQ(info);
    *info2 = -5; TaoFunctionReturn(0);
  }
  else if (neP->stepmax < neP->stepmin) {
    info = PetscInfo2(tao,"TaoApply_BoundLineSearch:Line search error: stepmax (%g) < stepmin (%g)\n", neP->stepmax,neP->stepmin); CHKERRQ(info);
    *info2 = -6; TaoFunctionReturn(0);
  } 
  else if (neP->maxfev < 0.0) {
    info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Line search error: maxfev (%d) < 0\n",neP->maxfev); CHKERRQ(info);
    *info2 = -7; TaoFunctionReturn(0);
  }

  /* Compute step length needed to make all variables equal a bound */ 
  /* Compute the smallest steplength that will make one nonbinding  */
  /* variable equal the bound */ 
  double unBoundStepNorm, boundStepNorm;
  info = TaoGetVariableBounds(tao,&XL,&XU); CHKERRQ(info);
  info = S->Norm2(&unBoundStepNorm); CHKERRQ(info);
  info = S->Negate(); CHKERRQ(info);
  info = S->BoundGradientProjection(S,XL,X,XU); CHKERRQ(info);
  info = S->Negate(); CHKERRQ(info);
  info = S->Norm2(&boundStepNorm); CHKERRQ(info);
  info = X->StepBoundInfo(XL,XU,S,&bstepmin1,&bstepmin2,&bstepmax); CHKERRQ(info);
  neP->stepmax=TaoMin(bstepmax,1.0e+15);

  info = PetscInfo2(tao,"TaoApply_BoundLineSearch:monitor: UnBoundNorm %22.12e, BoundNorm %22.12e \n",unBoundStepNorm,boundStepNorm); CHKERRQ(info);
  /* Check that search direction is a descent direction */

#if defined(PETSC_USE_COMPLEX)
  info = G->Dot(S,&cdginit);CHKERRQ(info); dginit = TaoReal(cdginit);
#else
  info = G->Dot(S,&dginit);CHKERRQ(info);
#endif

  if (dginit >= 0.0) {
    info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Search direction not a descent direction, dginit %22.12e \n",dginit); CHKERRQ(info);
    *info2 = 7; TaoFunctionReturn(0);
  }

  /* Initialization */
  neP->bracket = 0;
  *info2  = 0;
  stage1  = 1;
  finit   = *f;
  dgtest  = neP->ftol * dginit;
  width   = neP->stepmax - neP->stepmin;
  width1  = width * 2.0;
  info = W->CopyFrom(X);CHKERRQ(info);
  /* Variable dictionary:  
     stx, fx, dgx - the step, function, and derivative at the best step
     sty, fy, dgy - the step, function, and derivative at the other endpoint 
                   of the interval of uncertainty
     step, f, dg - the step, function, and derivative at the current step */

  stx = 0.0;
  fx  = finit;
  dgx = dginit;
  sty = 0.0;
  fy  = finit;
  dgy = dginit;
 
  neP->nfev = 0;
  tao->new_search=TAO_TRUE;
  for (i=0; i< neP->maxfev; i++) {
    /* Set min and max steps to correspond to the interval of uncertainty */
    if (neP->bracket) {
      neP->stepmin = TaoMin(stx,sty); 
      neP->stepmax = TaoMax(stx,sty); 
    } else {
      neP->stepmin = stx;
      neP->stepmax = *step + xtrapf * (*step - stx);
    }

    /* Force the step to be within the bounds */
    *step = TaoMax(*step,neP->stepmin);
    *step = TaoMin(*step,neP->stepmax);

    /* If an unusual termination is to occur, then let step be the lowest
       point obtained thus far */
    if (((neP->bracket) && (*step <= neP->stepmin || *step >= neP->stepmax)) ||
        ((neP->bracket) && (neP->stepmax - neP->stepmin <= neP->rtol * neP->stepmax)) ||
        (neP->nfev >= neP->maxfev - 1) || (neP->infoc == 0)) {
      *step = stx;
    }
 
#if defined(PETSC_USE_COMPLEX)
    cstep = *step;
    info = W->Waxpby(cstep,S,1.0,X);CHKERRQ(info);
#else
    info = W->Waxpby(*step,S,1.0,X);CHKERRQ(info); 	/* W = X + step*S */
#endif

    info=W->Median(XL,W,XU);CHKERRQ(info);
    tao->current_step=*step;
    if (tao->MeritFunctionGTSApply) {
      info = TaoComputeMeritFunctionGTS(tao,W,S,f,&dg); CHKERRQ(info);
      need_gradient = 1;
    } else {
      info = TaoComputeMeritFunctionGradient(tao,W,f,G);CHKERRQ(info);
      tao->new_search=TAO_FALSE;
      neP->nfev++;
#if defined(PETSC_USE_COMPLEX)
      info = G->Dot(S,&cdg);CHKERRQ(info); dg = TaoReal(cdg);
#else
      info = G->Dot(X,&dg1);CHKERRQ(info);	        /* dg = G^T S */
      info = G->Dot(W,&dg2);CHKERRQ(info);	        /* dg = G^T S */
      dg = (dg2-dg1) / (*step);
#endif
    }
    if (0 == i) {
      *f_full = *f;
    }

    if ((*f != *f) || (dg1 != dg1) || (dg2 != dg2) || (dg != dg)) {
      // User provided compute function generated Not-a-Number, assume
      // domain violation and set function value and directional
      // derivative to infinity.
      *f = TAO_INFINITY;
      dg = TAO_INFINITY;
      dg1 = TAO_INFINITY;
      dg2 = TAO_INFINITY;
    }

    ftest1 = finit + (*step) * dgtest;
    ftest2 = finit + (*step) * dgtest * neP->ftol;	// Armijo

    info = PetscInfo7(tao, "TaoApply_BoundLineSearch:monitor: function %22.15e ftest1 %22.15e ftest2 %22.15e step  %22.15e dg %22.15e gtol*dginit %22.15e bstepmin2 %22.15e \n",*f,ftest1,ftest2,*step,TaoAbsDouble(dg),neP->gtol*(-dginit),bstepmin2); CHKERRQ(info);

    /* Convergence testing */
    if ((*f <= ftest1) && (TaoAbsDouble(dg) <= neP->gtol*(-dginit))) {
      info = PetscInfo(tao, "TaoApply_BoundLineSearch:Line search success: Sufficient decrease and directional deriv conditions hold\n"); CHKERRQ(info);
      *info2 = 0; 
      break;
    }

    /* Check Armijo if beyond the first breakpoint */
    if ((*f <= ftest2) && (*step >= bstepmin2)) {
      info = PetscInfo(tao,"TaoApply_BoundLineSearch:Line search success: Sufficient decrease\n"); CHKERRQ(info);
      *info2 = 0; 
      break;
    }

    /* Checks for bad cases */
    if (((neP->bracket) && (*step <= neP->stepmin||*step >= neP->stepmax)) || (!neP->infoc)) {
      info = PetscInfo(tao,"TaoApply_LineSearch:Rounding errors may prevent further progress.  May not be a step satisfying\n"); CHKERRQ(info);
      info = PetscInfo(tao,"TaoApply_BoundLineSearch:sufficient decrease and curvature conditions. Tolerances may be too small.\n"); CHKERRQ(info);
      *info2 = 6; break;
    }
    if ((*step == neP->stepmax) && (*f <= ftest1) && (dg <= dgtest)) {
      info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Step is at the upper bound, stepmax (%g)\n",neP->stepmax); CHKERRQ(info);
      *info2 = 5; break;
    }
    if ((*step == neP->stepmin) && (*f >= ftest1) && (dg >= dgtest)) {
      info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Step is at the lower bound, stepmin (%g)\n",neP->stepmin); CHKERRQ(info);
      *info2 = 4; break;
    }
    if (neP->nfev >= neP->maxfev) {
      info = PetscInfo2(tao,"TaoApply_BoundLineSearch:Number of line search function evals (%d) > maximum (%d)\n",neP->nfev,neP->maxfev); CHKERRQ(info);
      *info2 = 3; break;
    }
    if ((neP->bracket) && (neP->stepmax - neP->stepmin <= neP->rtol*neP->stepmax)){
      info = PetscInfo1(tao,"TaoApply_BoundLineSearch:Relative width of interval of uncertainty is at most rtol (%g)\n",neP->rtol); CHKERRQ(info);
      *info2 = 2; break;
    }

    /* In the first stage, we seek a step for which the modified function
        has a nonpositive value and nonnegative derivative */
    if ((stage1) && (*f <= ftest1) && (dg >= dginit * TaoMin(neP->ftol, neP->gtol)))
      stage1 = 0;

    /* A modified function is used to predict the step only if we
       have not obtained a step for which the modified function has a 
       nonpositive function value and nonnegative derivative, and if a
       lower function value has been obtained but the decrease is not
       sufficient */

    if ((stage1) && (*f <= fx) && (*f > ftest1)) {
      fm   = *f - *step * dgtest;	/* Define modified function */
      fxm  = fx - stx * dgtest;	        /* and derivatives */
      fym  = fy - sty * dgtest;
      dgm  = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

      /* Update the interval of uncertainty and compute the new step */
      info = TaoStep_LineSearch(tao,&stx,&fxm,&dgxm,&sty,&fym,&dgym,step,&fm,&dgm);CHKERRQ(info);

      fx  = fxm + stx * dgtest;	/* Reset the function and */
      fy  = fym + sty * dgtest;	/* gradient values */
      dgx = dgxm + dgtest; 
      dgy = dgym + dgtest; 
    } else {
      /* Update the interval of uncertainty and compute the new step */
      info = TaoStep_LineSearch(tao,&stx,&fx,&dgx,&sty,&fy,&dgy,step,f,&dg);CHKERRQ(info);
    }

    /* Force a sufficient decrease in the interval of uncertainty */
    if (neP->bracket) {
      if (TaoAbsDouble(sty - stx) >= 0.66 * width1) *step = stx + 0.5*(sty - stx);
      width1 = width;
      width = TaoAbsDouble(sty - stx);
    }
  }

  /* Finish computations */
  info = X->CopyFrom(W); CHKERRQ(info);
  if (need_gradient) {
    info = TaoComputeMeritGradient(tao,X,G); CHKERRQ(info);
  }
  info = PetscInfo2(tao,"TaoApply_BoundLineSearch:%d function evals in line search, step = %10.4e\n",neP->nfev,*step); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_LineSearch"
static int TaoDestroy_LineSearch(TAO_SOLVER tao, void *ctx)
{
  int  info;

  TaoFunctionBegin;
  info = TaoFree(tao->linectx);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_LineSearch"
static int TaoSetOptions_LineSearch(TAO_SOLVER tao, void*linectx)
{
  TAO_LINESEARCH *ctx = (TAO_LINESEARCH *)linectx;
  int            info;

  TaoFunctionBegin;
  info = TaoOptionsHead("More-Thuente line line search options for minimization");CHKERRQ(info);
  info = TaoOptionInt("-tao_ls_maxfev","max function evals in line search","",ctx->maxfev,&ctx->maxfev,0);CHKERRQ(info);
  info = TaoOptionDouble("-tao_ls_ftol","tol for sufficient decrease","",ctx->ftol,&ctx->ftol,0);CHKERRQ(info);
  info = TaoOptionDouble("-tao_ls_gtol","tol for curvature condition","",ctx->gtol,&ctx->gtol,0);CHKERRQ(info);
  info = TaoOptionDouble("-tao_ls_rtol","relative tol for acceptable step","",ctx->rtol,&ctx->rtol,0);CHKERRQ(info);
  info = TaoOptionDouble("-tao_ls_stepmin","lower bound for step","",ctx->stepmin,&ctx->stepmin,0);CHKERRQ(info);
  info = TaoOptionDouble("-tao_ls_stepmax","upper bound for step","",ctx->stepmax,&ctx->stepmax,0);CHKERRQ(info);
  info = TaoOptionsTail();CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoView_LineSearch"
static int TaoView_LineSearch(TAO_SOLVER tao, void*ctx)
{
  TAO_LINESEARCH *ls = (TAO_LINESEARCH *)ctx;
  int            info;

  TaoFunctionBegin;
  info=TaoPrintInt(tao,"    More'-Thuente Line search: maxf=%d,",ls->maxfev);CHKERRQ(info);
  info=TaoPrintDouble(tao," ftol=%g,",ls->ftol);CHKERRQ(info);
  info=TaoPrintDouble(tao," rtol=%g,",ls->rtol);CHKERRQ(info);
  info=TaoPrintDouble(tao," gtol=%g\n",ls->gtol);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCreateMoreThuenteLineSearch"
/*@C
   TaoCreateMoreThuenteLineSearch - Create a line search

   Input Parameters:
+  tao - TAO_SOLVER context
.  fftol - the sufficient descent parameter , greater than 0.
-  ggtol - the curvature tolerance, greater than 0, less than 1.


   Note:
   If either fftol or ggtol is 0, default parameters will be used.

   Note:
   This algorithm is taken from More' and Thuente, "Line search algorithms
   with guaranteed sufficient decrease", Argonne National Laboratory, 
   Technical Report MCS-P330-1092.

   Note:
   This line search enforces the strong Wolfe conditions for unconstrained
   optimization.  This routine is used within the following TAO unconstrained
   minimization solvers: Newton linesearch (tao_nls), limited memory variable 
   metric (tao_lmvm), and nonlinear conjugate gradient methods.

   Level: developer

.keywords: TAO_SOLVER, linesearch
@*/
int TaoCreateMoreThuenteLineSearch(TAO_SOLVER tao, double fftol, double ggtol)
{
  int info;
  TAO_LINESEARCH *neP;

  TaoFunctionBegin;

  info = TaoNew(TAO_LINESEARCH,&neP); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_LINESEARCH)); CHKERRQ(info);

  if (fftol<=0) fftol=0.0001;
  if (ggtol<=0) ggtol=0.9;

  neP->ftol		  = fftol;
  neP->rtol		  = 1.0e-10;
  neP->gtol		  = ggtol;
  neP->stepmin		  = 1.0e-20;
  neP->stepmax		  = 1.0e+20;
  neP->nfev		  = 0; 
  neP->bracket		  = 0; 
  neP->infoc              = 1;
  neP->maxfev		  = 30;

  info = TaoSetLineSearch(tao,0,
			  TaoSetOptions_LineSearch,
			  TaoApply_LineSearch,
			  TaoView_LineSearch,
			  TaoDestroy_LineSearch,
			  (void *) neP);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCreateMoreThuenteBoundLineSearch"
/*@C
   TaoCreateMoreThuenteBoundLineSearch - Create a line search

   Input Parameters:
+  tao - TAO_SOLVER context
.  fftol - the sufficient descent parameter , greater than 0.
-  ggtol - the curvature tolerance, greater than 0, less than 1.


   Note:
   If either fftol or ggtol is 0, default parameters will be used.

   Note:
   This algorithm is is a modification of the algorithm by More' and Thuente.
   The modifications concern bounds.  This algorithm steps in the direction
   passed into this routine.  This point get projected back into the feasible set.
   It tries to satisfy the Wolfe conditions, but in the context of bound constrained
   optimization, there may not be a point in the piecewise linear 
   path that satisfies the Wolfe conditions.  When the active set
   is changing, decrease in the objective function may be sufficient
   to terminate this line search.

   Note:
   This routine is used within the following TAO bound constrained
   minimization solvers: Newton trust region (tao_tron) and limited memory variable 
   metric (tao_blmvm).

   Level: developer

.keywords: TAO_SOLVER, linesearch
@*/
int TaoCreateMoreThuenteBoundLineSearch(TAO_SOLVER tao, double fftol, double ggtol)
{
  int info;
  TAO_LINESEARCH *neP;

  TaoFunctionBegin;

  info = TaoNew(TAO_LINESEARCH,&neP); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_LINESEARCH)); CHKERRQ(info);

  if (fftol<=0) fftol=0.0001;
  if (ggtol<=0) ggtol=0.9;

  neP->ftol		  = fftol;
  neP->rtol		  = 1.0e-10;
  neP->gtol		  = ggtol;
  neP->stepmin		  = 1.0e-20;
  neP->stepmax		  = 1.0e+20;
  neP->nfev		  = 0; 
  neP->bracket		  = 0; 
  neP->infoc              = 1;
  neP->maxfev		  = 30;

  info = TaoSetLineSearch(tao,0,
			  TaoSetOptions_LineSearch,
			  TaoApply_BoundLineSearch,
			  TaoView_LineSearch,
			  TaoDestroy_LineSearch,
			  (void *) neP);CHKERRQ(info);

  TaoFunctionReturn(0);
}

