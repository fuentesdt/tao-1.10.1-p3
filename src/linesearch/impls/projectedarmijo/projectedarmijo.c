#include "projectedarmijo.h"

#define REPLACE_FIFO 1
#define REPLACE_MRU  2

#define REFERENCE_MAX 1
#define REFERENCE_AVE 2
#define REFERENCE_MEAN 3

#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_ProjectedArmijo"
static int TaoDestroy_ProjectedArmijo(TAO_SOLVER tao, void*ctx)
{
  TAO_PROJECTEDARMIJO *ls = (TAO_PROJECTEDARMIJO *)ctx;
  int info;

  TaoFunctionBegin;
  if (ls->work != TAO_NULL) {
    delete ls->work;
  }

  if (ls->memory != TAO_NULL) {
    info = TaoFree(ls->memory); CHKERRQ(info);
    ls->memory = TAO_NULL;
  }
  info = TaoFree(ls); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_ProjectedArmijo"
static int TaoSetOptions_ProjectedArmijo(TAO_SOLVER tao, void*ctx)
{
  TAO_PROJECTEDARMIJO *ls = (TAO_PROJECTEDARMIJO *)ctx;
  int info;

  TaoFunctionBegin;
  info = TaoOptionsHead("Projected Armijo linesearch options");CHKERRQ(info);
  info = TaoOptionDouble("-tao_projected_armijo_alpha", "initial reference constant", "", ls->alpha, &ls->alpha, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_projected_armijo_beta", "decrease constant", "", ls->beta, &ls->beta, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_projected_armijo_sigma", "acceptance constant", "", ls->sigma, &ls->sigma, 0); CHKERRQ(info);
  info = TaoOptionInt("-tao_projected_armijo_memory_size", "number of historical elements", "", ls->memorySize, &ls->memorySize, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_projected_armijo_minimum_step", "minimum acceptable step", "", ls->minimumStep, &ls->minimumStep, 0); CHKERRQ(info);
  info = TaoOptionInt("-tao_projected_armijo_reference_policy", "policy for updating reference value", "", ls->referencePolicy, &ls->referencePolicy, 0); CHKERRQ(info);
  info = TaoOptionInt("-tao_projected_armijo_replacement_policy", "policy for updating memory", "", ls->replacementPolicy, &ls->replacementPolicy, 0); CHKERRQ(info);
  info = TaoOptionsTail();CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoView_ProjectedArmijo"
static int TaoView_ProjectedArmijo(TAO_SOLVER tao, void *ctx)
{
  TAO_PROJECTEDARMIJO *ls = (TAO_PROJECTEDARMIJO *)ctx;
  int info;

  TaoFunctionBegin;

  info=TaoPrintDouble(tao,"  Projected Armijo linesearch: alpha=%g",ls->alpha);CHKERRQ(info);
  info=TaoPrintDouble(tao," beta=%g ",ls->beta);CHKERRQ(info);
  info=TaoPrintDouble(tao,"sigma=%g ",ls->sigma);CHKERRQ(info);
  info=TaoPrintDouble(tao,"minstep=%g,",ls->minimumStep);CHKERRQ(info);
  info=TaoPrintInt(tao,"memsize=%d\n",ls->memorySize);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoApply_PreProjectedArmijo"
static int TaoApply_PreProjectedArmijo(TAO_SOLVER tao, TAO_PROJECTEDARMIJO *ls,
				       double f, double step,
				       double *ref, TaoInt *idx, TaoInt *info2)
{
  int info;
  TaoInt i;

  TaoFunctionBegin;

  *info2 = 0;

  // Check linesearch parameters
  if (step < 0) {
    info = PetscInfo1(tao, "TaoApply_ProjectedArmijo:Line search error: step (%g) < 0\n", step); CHKERRQ(info);
    *info2 = -1; 
    TaoFunctionReturn(0);
  } else if (ls->alpha < 1) {
    info = PetscInfo1(tao,"TaoApply_ProjectedArmijo:Line search error: alpha (%g) < 1\n", ls->alpha); CHKERRQ(info);
    *info2 = -2; 
    TaoFunctionReturn(0);
  } else if ((ls->beta <= 0) || (ls->beta >= 1)) {
    info = PetscInfo1(tao,"TaoApply_ProjectedArmijo:Line search error: beta (%g) invalid\n", ls->beta); CHKERRQ(info);
    *info2 = -3; 
    TaoFunctionReturn(0);
  } else if ((ls->sigma <= 0) || (ls->sigma >= 0.5)) {
    info = PetscInfo1(tao,"TaoApply_ProjectedArmijo:Line search error: sigma (%g) invalid\n", ls->sigma); CHKERRQ(info);
    *info2 = -4; 
    TaoFunctionReturn(0);
  } else if (ls->minimumStep <= 0) {
    info = PetscInfo1(tao,"TaoApply_ProjectedArmijo:Line search error: minimum_step (%g) <= 0\n", ls->minimumStep); CHKERRQ(info);
    *info2 = -5; 
    TaoFunctionReturn(0);
  } else if (ls->memorySize < 1) {
    info = PetscInfo1(tao,"TaoApply_ProjectedArmijo:Line search error: memory_size (%d) < 1\n", ls->memorySize); CHKERRQ(info);
    *info2 = -6; 
    TaoFunctionReturn(0);
  } else if ((ls->referencePolicy != REFERENCE_MAX) &&
             (ls->referencePolicy != REFERENCE_AVE) &&
	     (ls->referencePolicy != REFERENCE_MEAN)){
    info = PetscInfo(tao,"TaoApply_ProjectedArmijo:Line search error: reference_policy invalid\n"); CHKERRQ(info);
    *info2 = -7; 
    TaoFunctionReturn(0);
  } else if ((ls->replacementPolicy != REPLACE_FIFO) && 
             (ls->replacementPolicy != REPLACE_MRU)) {
    info = PetscInfo(tao,"TaoApply_ProjectedArmijo:Line search error: replacement_policy invalid\n"); CHKERRQ(info);
    *info2 = -8; 
    TaoFunctionReturn(0);
  }

  // Check to see of the memory has been allocated.  If not, allocate
  // the historical array and populate it with the initial function
  // values.

  if (ls->memory == TAO_NULL) {
    info = TaoMalloc(sizeof(double)*ls->memorySize, &ls->memory);CHKERRQ(info);
    info = PetscLogObjectMemory(tao, sizeof(double)*ls->memorySize); CHKERRQ(info);
  }

  if (tao->iter == 0) {
    for (i = 0; i < ls->memorySize; i++) {
      ls->memory[i] = ls->alpha*(f);
    }

    ls->current = 0;
    ls->lastReference = ls->memory[0];
  }

  // Calculate reference value (MAX)
  *ref = ls->memory[0];
  *idx = 0;

  for (i = 1; i < ls->memorySize; i++) {
    if (ls->memory[i] > *ref) {
      *ref = ls->memory[i];
      *idx = i;
    }
  }

  if (ls->referencePolicy == REFERENCE_AVE) {
    *ref = 0;
    for (i = 0; i < ls->memorySize; i++) {
      *ref += ls->memory[i];
    }
    *ref = *ref / ls->memorySize;
    *ref = TaoMax(*ref, ls->memory[ls->current]);
  } else if (ls->referencePolicy == REFERENCE_MEAN) {
    *ref = TaoMin(*ref, 0.5*(ls->lastReference + ls->memory[ls->current]));
  }

  TaoFunctionReturn(0);
}

#undef __FUNCT__ 
#define __FUNCT__ "TaoApply_PostProjectedArmijo"
static int TaoApply_PostProjectedArmijo(TAO_SOLVER tao, TAO_PROJECTEDARMIJO *ls,
					double f, double step,
					double ref, TaoInt idx, TaoInt *info2)
{
  int info;
  TaoFunctionBegin;

  *info2 = 0;

  // Check termination
  if (step < ls->minimumStep) {
    info = PetscInfo(tao, "TaoApply_ProjectedArmijo:Step is at lower bound.\n"); CHKERRQ(info);
    *info2 = 1;
    TaoFunctionReturn(0);
  }

  // Successful termination, update memory
  ls->lastReference = ref;
  if (ls->replacementPolicy == REPLACE_FIFO) {
    ls->memory[ls->current++] = f;
    if (ls->current >= ls->memorySize) {
      ls->current = 0;
    }
  } else {
    ls->current = idx;
    ls->memory[idx] = f;
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_ProjectedArmijo"
/* @ TaoApply_ProjectedArmijo - This routine performs a linesearch. It
   backtracks until the (nonmonotone) Projected Armijo conditions are satisfied.

   Input Parameters:
+  tao - TAO_SOLVER context
.  X - current iterate (on output X contains new iterate, X + step*S)
.  S - search direction
.  f - merit function evaluated at X
.  G - gradient of merit function evaluated at X
.  W - work vector
-  step - initial estimate of step length

   Output parameters:
+  f - merit function evaluated at new iterate, X + step*S
.  G - gradient of merit function evaluated at new iterate, X + step*S
.  X - new iterate
-  step - final step length

   Info is set to one of:
.   0 - the line search succeeds; the sufficient decrease
   condition and the directional derivative condition hold

   negative number if an input parameter is invalid
-   -1 -  step < 0 

   positive number > 1 if the line search otherwise terminates
+    1 -  Step is at the lower bound, stepmin.
@ */

static int TaoApply_ProjectedArmijo(TAO_SOLVER tao, TaoVec *X, TaoVec *G, 
				    TaoVec *S, TaoVec *W, 
				    double *f, double *f_full, double *step,
				    TaoInt *info2, void *ctx)
{
  TAO_PROJECTEDARMIJO *ls = (TAO_PROJECTEDARMIJO *)ctx;
  TaoVec *L, *U, *work;
  double ref, innerd, t;
  TaoInt idx;
  int  info;
  TaoTruth flag;

  TaoFunctionBegin;

  info = TaoApply_PreProjectedArmijo(tao, ls, *f, *step, &ref, &idx, info2);
  if (*info2) {
    TaoFunctionReturn(0);
  }

  if (ls->work!=TAO_NULL){
    info=X->Compatible(ls->work,&flag); CHKERRQ(info);
    if (flag==TAO_FALSE){
      info=TaoVecDestroy(ls->work); CHKERRQ(info);
      ls->work=TAO_NULL;
    }
  }

  if (ls->work == TAO_NULL) {
     G->Clone(&(ls->work));
  }

  info = TaoGetVariableBounds(tao, &L, &U);
  work = ls->work;

  const double sigma = ls->sigma;
  const double beta = ls->beta;

  t = *step;
  tao->new_search=TAO_TRUE;
  while (t >= ls->minimumStep) {
    // Calculate iterate
    info = W->Waxpby(1.0, X, t, S); CHKERRQ(info);
    info = W->PointwiseMaximum(W, L); CHKERRQ(info);
    info = W->PointwiseMinimum(W, U); CHKERRQ(info);

    info = work->Waxpby(1.0, X, -1.0, W); CHKERRQ(info);
    info = work->Dot(G, &innerd); CHKERRQ(info);

    if (innerd > 0) {
      // Calculate function at new iterate
      tao->current_step=t;
      info = TaoComputeMeritFunction(tao, W, f); CHKERRQ(info);
      tao->new_search=TAO_FALSE;
      if (*step == t) {
        *f_full = *f;
      }

      // Check descent condition
      if (*f <= ref - sigma*innerd) {
        break;
      }
    }
    else if (*step == t) {
      tao->current_step=t;
      info = TaoComputeMeritFunction(tao, W, f_full); CHKERRQ(info);
      tao->new_search=TAO_FALSE;
    }

    t *= beta;
  }

  info = TaoApply_PostProjectedArmijo(tao, ls, *f, t, ref, idx, info2);

  // Update iterate and compute gradient
  *step = t;
  info = X->CopyFrom(W); CHKERRQ(info);
  tao->current_step=t;
  info = TaoComputeMeritGradient(tao, X, G); CHKERRQ(info);

  // Finish computations
  info = PetscInfo1(tao,"TaoApply_ProjectedArmijo:step = %10.4f\n",*step); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_NDProjectedArmijo"
/* @ TaoApply_NDProjectedArmijo - This routine performs a linesearch. It
   backtracks until the (nonmonotone) Projected Armijo conditions are 
   satisfied.  This is a modified version for a nondifferentiable function.

   Input Parameters:
+  tao - TAO_SOLVER context
.  X - current iterate (on output X contains new iterate, X + step*S)
.  S - search direction
.  f - merit function evaluated at X
-  step - initial estimate of step length

   Output parameters:
+  f - merit function evaluated at new iterate, X + step*S
.  X - new iterate
-  step - final step length

   Info is set to one of:
.   0 - the line search succeeds; the sufficient decrease
   condition and the directional derivative condition hold

   negative number if an input parameter is invalid
-   -1 -  step < 0 

   positive number > 1 if the line search otherwise terminates
+    1 -  Step is at the lower bound, stepmin.
@ */

static int TaoApply_NDProjectedArmijo(TAO_SOLVER tao, TaoVec *X, TaoVec *G, 
				      TaoVec *S, TaoVec *W, 
				      double *f, double *f_full, double *step,
				      TaoInt *info2, void *ctx)
{
  TAO_PROJECTEDARMIJO *ls = (TAO_PROJECTEDARMIJO *)ctx;
  TaoVec *L, *U;
  double ref, t;
  int info;
  TaoInt idx;

  TaoFunctionBegin;

  info = TaoApply_PreProjectedArmijo(tao, ls, *f, *step, &ref, &idx, info2);
  if (*info2) {
    TaoFunctionReturn(0);
  }

  info = TaoGetVariableBounds(tao, &L, &U);

  const double sigma = ls->sigma;
  const double beta = ls->beta;

  t = *step;
  tao->new_search=TAO_TRUE;
  while (t >= ls->minimumStep) {
    // Calculate iterate
    info = W->Waxpby(1.0, X, t, S); CHKERRQ(info);
    info = W->PointwiseMaximum(W, L); CHKERRQ(info);
    info = W->PointwiseMinimum(W, U); CHKERRQ(info);

    // Calculate function at new iterate

    tao->current_step=t;
    info = TaoComputeMeritFunction(tao, W, f); CHKERRQ(info);
    tao->new_search=TAO_FALSE;
    if (*step == t) {
      *f_full = *f;
    }

    // Check descent condition
    if (*f <= (1 - sigma*t)*ref) {
        break;
    }
    
    t *= beta;
  }

  info = TaoApply_PostProjectedArmijo(tao, ls, *f, t, ref, idx, info2);

  // Update iterate and compute gradient
  *step = t;
  info = X->CopyFrom(W); CHKERRQ(info);
  tao->current_step=t;
  info = TaoComputeMeritGradient(tao, X, G); CHKERRQ(info);


  // Finish computations
  info = PetscInfo1(tao,"TaoApply_NDProjectedArmijo:step = %10.4f\n",*step); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCreateProjectedArmijoLineSearch"
/*@
   TaoCreateProjectedArmijoLineSearch - Create a non-monotone projected linesearch

   Input Parameters:
.  tao - TAO_SOLVER context


   Note:
   This algorithm is taken from the following references -- 

   Armijo, "Minimization of Functions Having Lipschitz Continuous
     First-Partial Derivatives," Pacific Journal of Mathematics, volume 16,
     pages 1-3, 1966.
   Ferris and Lucidi, "Nonmonotone Stabilization Methods for Nonlinear
     Equations," Journal of Optimization Theory and Applications, volume 81,
     pages 53-71, 1994.
   Grippo, Lampariello, and Lucidi, "A Nonmonotone Line Search Technique
     for Newton's Method," SIAM Journal on Numerical Analysis, volume 23,
     pages 707-716, 1986.
   Grippo, Lampariello, and Lucidi, "A Class of Nonmonotone Stabilization
     Methods in Unconstrained Optimization," Numerische Mathematik, volume 59,
     pages 779-805, 1991.

   Note:
   This line seach enforces non-monotone Armijo descent conditions for
   bounds constrained optimization.  This routine is used within the 
   following TAO solvers: feasible semismooth with linesearch (tao_ssfls).

   Level: developer

.keywords: TAO_SOLVER, linesearch
@*/
int TaoCreateProjectedArmijoLineSearch(TAO_SOLVER tao)
{
  TAO_PROJECTEDARMIJO *ls;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_PROJECTEDARMIJO,&ls); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_PROJECTEDARMIJO)); CHKERRQ(info);

  ls->work = TAO_NULL;
  ls->memory = TAO_NULL;
  ls->alpha = 1.0;
  ls->beta = 0.5;
  ls->sigma = 1e-4;
  ls->minimumStep = TAO_EPSILON;
  ls->memorySize = 1;
  ls->referencePolicy = REFERENCE_MAX;
  ls->replacementPolicy = REPLACE_MRU;

  info = TaoSetLineSearch(tao,0,
			  TaoSetOptions_ProjectedArmijo,
			  TaoApply_ProjectedArmijo,
			  TaoView_ProjectedArmijo,
			  TaoDestroy_ProjectedArmijo,
			  (void *) ls);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCreateNDProjectedArmijoLineSearch"
/*@
   TaoCreateNDProjectedArmijoLineSearch - Create a non-monotone projected linesearch
     for a nondifferentiable function

   Input Parameters:
.  tao - TAO_SOLVER context


   Note:
   This algorithm is taken from the following references -- 

   Armijo, "Minimization of Functions Having Lipschitz Continuous
     First-Partial Derivatives," Pacific Journal of Mathematics, volume 16,
     pages 1-3, 1966.
   Ferris and Lucidi, "Nonmonotone Stabilization Methods for Nonlinear
     Equations," Journal of Optimization Theory and Applications, volume 81,
     pages 53-71, 1994.
   Grippo, Lampariello, and Lucidi, "A Nonmonotone Line Search Technique
     for Newton's Method," SIAM Journal on Numerical Analysis, volume 23,
     pages 707-716, 1986.
   Grippo, Lampariello, and Lucidi, "A Class of Nonmonotone Stabilization
     Methods in Unconstrained Optimization," Numerische Mathematik, volume 59,
     pages 779-805, 1991.

   Note:
   This line seach enforces non-monotone Armijo descent conditions for
   bounds constrained optimization.  This routine is used within the 
   following TAO solvers: feasible semismooth with linesearch (tao_ssfls).

   Level: developer

.keywords: TAO_SOLVER, linesearch
@*/
int TaoCreateNDProjectedArmijoLineSearch(TAO_SOLVER tao)
{
  TAO_PROJECTEDARMIJO *ls;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_PROJECTEDARMIJO,&ls); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_PROJECTEDARMIJO));CHKERRQ(info);

  ls->work = TAO_NULL;
  ls->memory = TAO_NULL;
  ls->alpha = 1.0;
  ls->beta = 0.5;
  ls->sigma = 1e-4;
  ls->minimumStep = TAO_EPSILON;
  ls->memorySize = 1;
  ls->referencePolicy = REFERENCE_MAX;
  ls->replacementPolicy = REPLACE_MRU;

  info = TaoSetLineSearch(tao,0,
			  TaoSetOptions_ProjectedArmijo,
			  TaoApply_NDProjectedArmijo,
			  TaoView_ProjectedArmijo,
			  TaoDestroy_ProjectedArmijo,
			  (void *) ls);CHKERRQ(info);

  TaoFunctionReturn(0);
}

