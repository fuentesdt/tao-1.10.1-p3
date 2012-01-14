
#include "src/bound/impls/gpcg/gpcglinesearch.h"    /*I "tao_solver.h" I*/

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_LineSearch2"
static int TaoGPCGDestroyLineSearch(TAO_SOLVER tao, void*lsctx)
{
  int  info;
  TAO_GPCGLINESEARCH *ctx = (TAO_GPCGLINESEARCH *)lsctx;

  TaoFunctionBegin;
  if (ctx->setupcalled==1){
    info = TaoVecDestroy(ctx->W2);CHKERRQ(info);
    info = TaoVecDestroy(ctx->Gold);CHKERRQ(info);
  }
  info = TaoFree(ctx);CHKERRQ(info);
  TaoFunctionReturn(0);
}
/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_LineSearch2"
static int TaoGPCGSetOptionsLineSearch(TAO_SOLVER tao, void*linectx)
{
  TAO_GPCGLINESEARCH *ctx = (TAO_GPCGLINESEARCH *)linectx;
  double         tmp;
  int            info;
  TaoInt         itmp;
  TaoTruth     flg;

  TaoFunctionBegin;
  info = TaoOptionsHead("GPCG line search options");CHKERRQ(info);

  info = TaoOptionInt("-tao_nls_maxfev","max function evals in line search",0,ctx->maxfev,&itmp,&flg);CHKERRQ(info);
  if (flg) {ctx->maxfev = itmp;}
  info = TaoOptionDouble("-tao_nls_ftol","tol for sufficient decrease",0,ctx->ftol,&tmp,&flg);CHKERRQ(info);
  if (flg) {ctx->ftol = tmp;}
  info = TaoOptionDouble("-tao_nls_gtol","tol for curvature condition",0,ctx->gtol,&tmp,&flg);CHKERRQ(info);
  if (flg) {ctx->gtol = tmp;}
  info = TaoOptionDouble("-tao_nls_rtol","relative tol for acceptable step",0,ctx->rtol,&tmp,&flg);CHKERRQ(info);
  if (flg) {ctx->rtol = tmp;}
  info = TaoOptionDouble("-tao_nls_stepmin","lower bound for step",0,ctx->stepmin,&tmp,&flg);CHKERRQ(info);
  if (flg) {ctx->stepmin = tmp;}
  info = TaoOptionDouble("-tao_nls_stepmax","upper bound for step",0,ctx->stepmax,&tmp,&flg);CHKERRQ(info);
  if (flg) {ctx->stepmax = tmp;}
  info = TaoOptionsTail();CHKERRQ(info);

  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoGPCGViewLineSearch"
static int TaoGPCGViewLineSearch(TAO_SOLVER tao,void *ctx)
{
  TAO_GPCGLINESEARCH *ls = (TAO_GPCGLINESEARCH *)ctx;
  int            info;

  TaoFunctionBegin;
  info = TaoPrintInt(tao,"  Line search: maxf=%d,",ls->maxfev);CHKERRQ(info);
  info = TaoPrintDouble(tao," ftol=%g,",ls->ftol);CHKERRQ(info);
  info = TaoPrintDouble(tao," rtol=%g,",ls->rtol);CHKERRQ(info);
  info = TaoPrintDouble(tao," gtol=%g\n",ls->gtol);CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoGPCGApplyLineSearch"
static int TaoGPCGApplyLineSearch(TAO_SOLVER tao,TaoVec* X, 
			   TaoVec* G,TaoVec* S,TaoVec* W,
			   double *f, double *f_full, double *step, TaoInt *info2,
			   void*ctx)
{
  TAO_GPCGLINESEARCH *neP = (TAO_GPCGLINESEARCH *) ctx;
  int       info;
  TaoInt i;
  double zero=0.0;
  double d1,finit,actred,prered,rho, gdx;
  TaoVec* XL, *XU, *Xold=neP->W2,*Gold=neP->Gold;
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
  *info2=0;
  if (neP->setupcalled){
    info=X->Compatible(neP->W2,&flag); CHKERRQ(info);
    if (flag==TAO_FALSE){
      info=TaoVecDestroy(neP->W2); CHKERRQ(info);neP->W2=0;
      info=TaoVecDestroy(neP->Gold); CHKERRQ(info);neP->Gold=0;
      neP->setupcalled=0;
    }
  }

  if (neP->setupcalled==0){
    info = X->Clone(&neP->W2); CHKERRQ(info);
    Xold=neP->W2;
    info = X->Clone(&neP->Gold); CHKERRQ(info);
    Gold=neP->Gold;
    neP->setupcalled=1;
  }

  info = G->Dot(S,&gdx); CHKERRQ(info);
  info = Xold->CopyFrom(X); CHKERRQ(info);
  info = Gold->CopyFrom(G); CHKERRQ(info);
  info = TaoGetVariableBounds(tao,&XL,&XU); CHKERRQ(info);
  info = X->StepBoundInfo(XL,XU,S,&rho,&actred,&d1);CHKERRQ(info);
  rho=0; actred=0;
  *step = TaoMin(*step,d1);

  if (*step < zero) {
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Line search error: step (%g) < 0\n",*step); CHKERRQ(info);
    *info2 = -1; TaoFunctionReturn(0);
  } else if (neP->ftol < zero) {
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Line search error: ftol (%g) < 0\n",neP->ftol); CHKERRQ(info);
    *info2 = -2; TaoFunctionReturn(0);
  } else if (neP->rtol < zero) {
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Line search error: rtol (%g) < 0\n",neP->rtol); CHKERRQ(info);
    *info2 = -3; TaoFunctionReturn(0);
  } else if (neP->gtol < zero) {
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Line search error: gtol (%g) < 0\n",neP->gtol); CHKERRQ(info);
    *info2 = -4; TaoFunctionReturn(0);
  } else if (neP->stepmin < zero) {
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Line search error: stepmin (%g) < 0\n",neP->stepmin); CHKERRQ(info);
    *info2 = -5; TaoFunctionReturn(0);
  } else if (neP->stepmax < neP->stepmin) {
    info = PetscInfo2(tao,"TaoGPCGApplyLineSearch:Line search error: stepmax (%g) < stepmin (%g)\n",neP->stepmax,neP->stepmin); CHKERRQ(info);
    *info2 = -6; TaoFunctionReturn(0);
  } else if (neP->maxfev < zero) {
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Line search error: maxfev (%d) < 0\n",neP->maxfev); CHKERRQ(info);
    *info2 = -7; TaoFunctionReturn(0);
  }


  /* Check that search direction is a descent direction */
  /*
  info = VecDot(G,S,&dginit);CHKERRQ(info);  / * dginit = G^T S * /
  if (dginit >= zero) {
    info = PetscLogInfo((tao,"TaoGPCGApplyLineSearch:Search direction not a descent direction\n")); CHKERRQ(info);
    *info2 = 7; PetscFunctionReturn(0);
  }
  */
  /* Initialization */
  neP->nfev = 0;
  finit = *f;
  for (i=0; i< neP->maxfev; i++) {
    
    /* Force the step to be within the bounds */
    *step = TaoMax(*step,neP->stepmin);
    *step = TaoMin(*step,neP->stepmax);
    
    info = X->Waxpby(*step,S,1.0,Xold); CHKERRQ(info);
    info = X->Median(XL,X,XU); CHKERRQ(info);

    info = TaoGPCGComputeFunctionGradient(tao, X, f, G); CHKERRQ(info);
    if (0 == i) {
      *f_full = *f;
    }

    actred = *f - finit;
    info = W->Waxpby(-1.0,Xold,1.0,X); CHKERRQ(info);
    info = W->Dot(Gold,&prered); CHKERRQ(info);
    if (fabs(prered)<1.0e-100) prered=1.0e-12;
    rho = actred/prered;
    /* 
       If sufficient progress has been obtained, accept the
       point.  Otherwise, backtrack. 
    */

    if (rho > neP->ftol){
      break;
    } else{
      *step = (*step)/2;
    }
  }

  /* Convergence testing */
  
  if (*step <= neP->stepmin || *step >= neP->stepmax) {
    *info2 = 6;
    info = PetscInfo(tao,"TaoGPCGApplyLineSearch:Rounding errors may prevent further progress.  May not be a step satisfying\n"); CHKERRQ(info);
    info = PetscInfo(tao,"TaoGPCGApplyLineSearch:sufficient decrease and curvature conditions. Tolerances may be too small.\n"); CHKERRQ(info);
  }
  if (*step == neP->stepmax) {
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Step is at the upper bound, stepmax (%g)\n",neP->stepmax); CHKERRQ(info);
    *info2 = 5;
  }
  if (*step == neP->stepmin) {
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Step is at the lower bound, stepmin (%g)\n",neP->stepmin); CHKERRQ(info);
    *info2 = 4;
  }
  if (neP->nfev >= neP->maxfev) {
    info = PetscInfo2(tao,"TaoGPCGApplyLineSearch:Number of line search function evals (%d) > maximum (%d)\n",neP->nfev,neP->maxfev); CHKERRQ(info);
    *info2 = 3;
  }
  if ((neP->bracket) && (neP->stepmax - neP->stepmin <= neP->rtol*neP->stepmax)){
    info = PetscInfo1(tao,"TaoGPCGApplyLineSearch:Relative width of interval of uncertainty is at most rtol (%g)\n",neP->rtol); CHKERRQ(info);
    *info2 = 2;
  }
  /*
  if ((*f <= ftest1) && (PetscAbsDouble(dg) <= neP->gtol*(-dginit))) {
    info = PetscLogInfo((tao,"TaoGPCGApplyLineSearch:Line search success: Sufficient decrease and directional deriv conditions hold\n")); CHKERRQ(info);
    *info2 = 1;
  }
  */
  
  /* Finish computations */
  info = PetscInfo2(tao,"TaoGPCGApplyLineSearch:%d function evals in line search, step = %10.4f\n",neP->nfev,*step); CHKERRQ(info);

  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoGPCGCreateLineSearch"
int TaoGPCGCreateLineSearch(TAO_SOLVER tao)
{
  int info;
  TAO_GPCGLINESEARCH *neP;

  TaoFunctionBegin;

  info = TaoNew(TAO_GPCGLINESEARCH,&neP);CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_GPCGLINESEARCH)); CHKERRQ(info);
  neP->ftol		  = 0.05;
  neP->rtol		  = 0.0;
  neP->gtol		  = 0.0;
  neP->stepmin		  = 1.0e-20;
  neP->stepmax		  = 1.0e+20;
  neP->nfev		  = 0; 
  neP->bracket		  = 0; 
  neP->infoc              = 1;
  neP->maxfev		  = 30;
  neP->setupcalled        = 0;

  info = TaoSetLineSearch(tao,0,
			  TaoGPCGSetOptionsLineSearch,
			  TaoGPCGApplyLineSearch,
			  TaoGPCGViewLineSearch,
			  TaoGPCGDestroyLineSearch,
			  (void *) neP);CHKERRQ(info);

  TaoFunctionReturn(0);
}

