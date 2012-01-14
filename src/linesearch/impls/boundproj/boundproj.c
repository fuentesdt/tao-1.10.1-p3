
#include "boundproj.h"    /*I "tao_solver.h" I*/

#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_LineSearch"
static int TaoDestroy_LineSearch(TAO_SOLVER tao,void *linectx)
{
  int  info;
  TAO_LINESEARCH2 *ctx = (TAO_LINESEARCH2 *)linectx;

  TaoFunctionBegin;
  if (ctx->setupcalled==1){
    info=TaoVecDestroy(ctx->W2);CHKERRQ(info);
  }
  info = TaoFree(ctx);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_LineSearch"
static int TaoSetOptions_LineSearch(TAO_SOLVER tao,void *linectx)
{
  TAO_LINESEARCH2 *ctx = (TAO_LINESEARCH2 *)linectx;
  int            info;

  TaoFunctionBegin;
  info = TaoOptionsHead("Projected line search options");CHKERRQ(info);
  info = TaoOptionInt("-tao_ls_maxfev","max function evals in line search","",ctx->maxfev,&ctx->maxfev,0);CHKERRQ(info);
  info = TaoOptionDouble("-tao_ls_ftol","tol for sufficient decrease","",ctx->ftol,&ctx->ftol,0);CHKERRQ(info);
  info = TaoOptionDouble("-tao_ls_stepmin","lower bound for step","",ctx->stepmin,&ctx->stepmin,0);CHKERRQ(info);
  info = TaoOptionDouble("-tao_ls_stepmax","upper bound for step","",ctx->stepmax,&ctx->stepmax,0);CHKERRQ(info);
  info = TaoOptionsTail();CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoView_LineSearch"
static int TaoView_LineSearch(TAO_SOLVER tao,void *ctx)
{
  TAO_LINESEARCH2 *ls = (TAO_LINESEARCH2 *)ctx;
  int            info;

  TaoFunctionBegin;
  info = TaoPrintInt(tao,"  Line search: maxf=%d,",ls->maxfev);CHKERRQ(info);
  info = TaoPrintDouble(tao," ftol=%g\n",ls->ftol);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_LineSearch"
/* @ TaoApply_LineSearch - This routine performs a line search. It
   backtracks until the Armijo conditions are satisfied.

   Input Parameters:
+  tao - TAO_SOLVER context
.  X - current iterate (on output X contains new iterate, X + step*S)
.  S - search direction
.  f - objective function evaluated at X
.  G - gradient evaluated at X
.  Gold - work vector
.  gdx - inner product of gradient and the direction of the first linear manifold being searched
-  step - initial estimate of step length

   Output parameters:
+  f - objective function evaluated at new iterate, X + step*S
.  G - gradient evaluated at new iterate, X + step*S
.  X - new iterate
-  step - final step length

   Info is set to one of:
.   0 - the line search succeeds; the sufficient decrease
   condition and the directional derivative condition hold

   negative number if an input parameter is invalid
.   -1 -  step < 0 
.   -2 -  ftol < 0 
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
+    7 -  Search direction is not a descent direction.

   Notes:
   This routine is used within the TAO_NLS method.
@ */
static int TaoApply_LineSearch(TAO_SOLVER tao,TaoVec* X,TaoVec* G,TaoVec* S,TaoVec* Gold,double *f,double *f_full,
                        double *step,TaoInt *info2,void*ctx)
{
  TAO_LINESEARCH2 *neP = (TAO_LINESEARCH2 *) ctx;
  int       info;
  TaoInt i;
  double zero=0.0;
  double finit,gdx,dginit,actred,prered,rho,d1=0,d2=0,d3=0;
  TaoVec* XL, *XU, *Xtmp=neP->W2;
  TaoTruth flag;

  TaoFunctionBegin;
  /* neP->stepmin - lower bound for step */
  /* neP->stepmax - upper bound for step */
  /* neP->rtol 	  - relative tolerance for an acceptable step */
  /* neP->ftol 	  - tolerance for sufficient decrease condition */
  /* neP->gtol 	  - tolerance for curvature condition */
  /* neP->nfev 	  - number of function evaluations */
  /* neP->maxfev  - maximum number of function evaluations */

  /* Check input parameters for errors */

  *info2 = 0;
    /* After 2 failures, Check that search direction is a descent direction */
    /* This test is not really sufficient.  */
  if (neP->setupcalled){
    info=X->Compatible(neP->W2,&flag); CHKERRQ(info);
    if (flag==TAO_FALSE){
      info=TaoVecDestroy(neP->W2); CHKERRQ(info);neP->W2=0;
      neP->setupcalled=0;
    }
  }
  if (neP->setupcalled==0){
    info=X->Clone(&neP->W2); CHKERRQ(info);
    Xtmp=neP->W2;
    neP->setupcalled=1;
  }
  info = TaoGetVariableBounds(tao,&XL,&XU); CHKERRQ(info);
  info = Gold->ScaleCopyFrom(-1.0,S); CHKERRQ(info);
  info = Gold->BoundGradientProjection(S,XL,X,XU); CHKERRQ(info);
  info = Gold->Dot(G,&dginit);CHKERRQ(info);  /* dginit = G^T S */
  dginit*=-1;
  gdx=dginit;
  if (dginit >= zero) {
    info = G->CopyFrom(Gold); CHKERRQ(info);
    info = PetscInfo(tao,"TaoApply_LineSearch2:Search direction not a descent direction\n"); CHKERRQ(info);
    *info2 = 7; TaoFunctionReturn(0);
  }

  info = Xtmp->CopyFrom(X); CHKERRQ(info);
  info = Gold->CopyFrom(G); CHKERRQ(info);
  info=TaoGetVariableBounds(tao,&XL,&XU); CHKERRQ(info);
  if (*step < zero) {
    info = PetscInfo1(tao,"TaoApply_LineSearch:Line search error: step (%g) < 0\n",*step); CHKERRQ(info);
    *info2 = -1; TaoFunctionReturn(0);
  } else if (neP->ftol < zero) {
    info = PetscInfo1(tao,"TaoApply_LineSearch:Line search error: ftol (%g) < 0\n",neP->ftol); CHKERRQ(info);
    *info2 = -2; TaoFunctionReturn(0);
  } else if (neP->maxfev < zero) {
    info = PetscInfo1(tao,"TaoApply_LineSearch:Line search error: maxfev (%d) < 0\n",neP->maxfev); CHKERRQ(info);
    *info2 = -7; TaoFunctionReturn(0);
  }

  /* Initialization */
  neP->nfev = 0;
  finit = *f;
  tao->new_search=TAO_TRUE;
  for (i=0; i< neP->maxfev; i++) {
    
    /* Force the step to be within the bounds */
    *step = TaoMax(*step,neP->stepmin);
    *step = TaoMin(*step,neP->stepmax);
    
    if (0==i){
      info = Xtmp->Axpy(*step,S);CHKERRQ(info);
      info = Xtmp->Median(XL,Xtmp,XU);CHKERRQ(info);
      tao->current_step=*step;
      info = TaoComputeMeritFunctionGradient(tao,Xtmp,f,G);
      tao->new_search=TAO_FALSE;
      *f_full = *f;

      neP->nfev++;
      if (info==0 && *f==*f){ /* Successful function evaluation */
	actred = *f - finit;
	rho = actred/( (*step) * (-TaoAbsScalar(gdx)) );
	if (actred < 0 && rho > neP->ftol){
	  break;
	} else{
	  info=X->StepBoundInfo(XL,XU,S,&d3,&d2,&d1);CHKERRQ(info);
	  info = PetscInfo1(tao,"stepmax = %10f\n",d1); CHKERRQ(info);

	  *step = TaoMin(*step,d1);
	  info = Gold->Dot(Xtmp,&d1); CHKERRQ(info);
	  info = Gold->Dot(X,&prered); CHKERRQ(info);
	  prered=d1-prered;
	  rho = actred/(-TaoAbsScalar(prered));
	  if (actred < 0 && rho > neP->ftol){
	    break;
	  } else{
	        *step = (*step)/2;
	  }
	}
      } else { /* Function could not be evaluated at this point */
	*step = (*step)*0.7;
      }

    } else {
      info = Xtmp->Waxpby(*step,S,1.0,X);CHKERRQ(info);
      info = Xtmp->Median(XL,Xtmp,XU);CHKERRQ(info);
      
      info = G->Waxpby(-1,X,1.0,Xtmp);CHKERRQ(info);
      info = G->Dot(Gold,&prered); CHKERRQ(info);
      tao->current_step=*step;
      info = TaoComputeMeritFunctionGradient(tao,Xtmp,f,G); CHKERRQ(info);
      tao->new_search=TAO_FALSE;
      neP->nfev++;
      if (info==0 && *f==*f){ /* Successful function evaluation */
	actred = *f - finit;
	rho = actred/(-TaoAbsScalar(prered));
	/* 
	   If sufficient progress has been obtained, accept the
	   point.   Prered could be positive or negative.  
	   Otherwise, backtrack. 
	*/

	if (actred < 0 && rho > neP->ftol){
	  info = PetscInfo4(tao,"TaoApply_LineSearch: steplength: %g, actual reduction: %g (hopefully positive), Predicted reduction: %g, rho: %g\n",*step,-actred,prered,rho); CHKERRQ(info);
	  break;
	} else {
	  info = PetscInfo4(tao,"TaoApply_LineSearch: steplength: %g, actual reduction: %g (hopefully positive), Predicted reduction: %g, rho: %g\n",*step,-actred,prered,rho); CHKERRQ(info);
	  if (i<5){
	    *step = (*step)/2;
	  } else {
	    *step = (*step)/2;
	  }
	}
      } else { /* Function could not be evaluated at this point */
	  info = PetscInfo(tao,"TaoApply_LineSearch: Problem in function evaluation\n"); CHKERRQ(info);
	*step = (*step)*0.7;
      }
    }

  }
  /* Convergence testing */
  
  if (*step <= neP->stepmin || *step >= neP->stepmax) {
    *info2 = 6;
    info = PetscInfo(tao,"TaoApply_LineSearch:Rounding errors may prevent further progress.  May not be a step satisfying\n"); CHKERRQ(info);
    info = PetscInfo(tao,"TaoApply_LineSearch:sufficient decrease and curvature conditions. Tolerances may be too small.\n"); CHKERRQ(info);
  }
  if (*step == neP->stepmax) {
    info = PetscInfo1(tao,"TaoApply_LineSearch:Step is at the upper bound, stepmax (%g)\n",neP->stepmax); CHKERRQ(info);
    *info2 = 5;
  }
  if (*step == neP->stepmin) {
    info = PetscInfo1(tao,"TaoApply_LineSearch:Step is at the lower bound, stepmin (%g)\n",neP->stepmin); CHKERRQ(info);
    *info2 = 4;
  }
  if (neP->nfev >= neP->maxfev) {
    info = PetscInfo2(tao,"TaoApply_LineSearch:Number of line search function evals (%d) > maximum (%d)\n",neP->nfev,neP->maxfev); CHKERRQ(info);
    *info2 = 3;
  }
  if ((neP->bracket) && (neP->stepmax - neP->stepmin <= neP->rtol*neP->stepmax)){
    info = PetscInfo1(tao,"TaoApply_LineSearch:Relative width of interval of uncertainty is at most rtol (%g)\n",neP->rtol); CHKERRQ(info);
    *info2 = 2;
  }
  /*
  if ((*f <= ftest1) && (TaoAbsDouble(dg) <= neP->gtol*(-dginit))) {
    info = PetscLogInfo((tao,"TaoApply_LineSearch:Line search success: Sufficient decrease and directional deriv conditions hold\n")); CHKERRQ(info);
    *info2 = 1;
  }
  */
  
  /* Finish computations */
  info = X->CopyFrom(Xtmp); CHKERRQ(info);
  info = PetscInfo2(tao,"TaoApply_LineSearch:%d function evals in line search, step = %10.4f\n",neP->nfev,*step); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCreateProjectedLineSearch"
/*@C
   TaoCreateProjectedLineSearch - Create a line search

   Input Parameters:
.  tao - TAO_SOLVER context

   Note:
   This routine is used within the following TAO bound constrained 
   minimization solvers: TRON (tao_tron) and limited memory variable metric
   (tao_blmvm). This line search projects points onto the feasible, bound
   constrained region.  It only enforces the Armijo descent condition.

   Level: developer

.keywords: TAO_SOLVER, linesearch
@*/
int TaoCreateProjectedLineSearch(TAO_SOLVER tao)
{
  int info;
  TAO_LINESEARCH2 *neP;

  TaoFunctionBegin;

  info = TaoNew(TAO_LINESEARCH2,&neP); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_LINESEARCH2)); CHKERRQ(info);

  neP->ftol		  = 0.001;
  neP->rtol		  = 0.0;
  neP->gtol		  = 0.0;
  neP->stepmin		  = 1.0e-30;
  neP->stepmax		  = 1.0e+20;
  neP->nfev		  = 0; 
  neP->bracket		  = 0; 
  neP->infoc              = 1;
  neP->maxfev		  = 30;
  neP->W2                 = TAO_NULL;
  neP->setupcalled        = 0;

  info = TaoSetLineSearch(tao,0,
			  TaoSetOptions_LineSearch,
			  TaoApply_LineSearch,
			  TaoView_LineSearch,
			  TaoDestroy_LineSearch,
			  (void *) neP);CHKERRQ(info);

  TaoFunctionReturn(0);
}

