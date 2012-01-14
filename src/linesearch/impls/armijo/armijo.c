#include "armijo.h"

#define REPLACE_FIFO 1
#define REPLACE_MRU  2

#define REFERENCE_MAX  1
#define REFERENCE_AVE  2
#define REFERENCE_MEAN 3

#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_Armijo"
static int TaoDestroy_Armijo(TAO_SOLVER tao, void*ctx)
{
  TAO_ARMIJO *ls = (TAO_ARMIJO *)ctx;
  int info;

  TaoFunctionBegin;

  if (ls->memory != TAO_NULL) {
    info = TaoFree(ls->memory); CHKERRQ(info);
    ls->memory = TAO_NULL;
  }
  info = TaoFree(ls); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_Armijo"
static int TaoSetOptions_Armijo(TAO_SOLVER tao, void* ctx)
{
  TAO_ARMIJO *ls = (TAO_ARMIJO *)tao->linectx;
  int info;

  TaoFunctionBegin;
  info = TaoOptionsHead("Armijo linesearch options");CHKERRQ(info);
  info = TaoOptionDouble("-tao_armijo_alpha", "initial reference constant", "", ls->alpha, &ls->alpha, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_armijo_beta_inf", "decrease constant one", "", ls->beta_inf, &ls->beta_inf, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_armijo_beta", "decrease constant", "", ls->beta, &ls->beta, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_armijo_sigma", "acceptance constant", "", ls->sigma, &ls->sigma, 0); CHKERRQ(info);
  info = TaoOptionInt("-tao_armijo_memory_size", "number of historical elements", "", ls->memorySize, &ls->memorySize, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_armijo_minimum_step", "minimum acceptable step", "", ls->minimumStep, &ls->minimumStep, 0); CHKERRQ(info);
  info = TaoOptionInt("-tao_projected_armijo_reference_policy", "policy for updating reference value", "", ls->referencePolicy, &ls->referencePolicy, 0); CHKERRQ(info);
  info = TaoOptionInt("-tao_projected_armijo_replacement_policy", "policy for updating memory", "", ls->replacementPolicy, &ls->replacementPolicy, 0); CHKERRQ(info);
  info = TaoOptionsTail();CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoView_Armijo"
static int TaoView_Armijo(TAO_SOLVER tao, void*ctx)
{
  TAO_ARMIJO *ls = (TAO_ARMIJO *)ctx;
  int info;

  TaoFunctionBegin;
  
  info=TaoPrintDouble(tao,"  Armijo linesearch: alpha=%g",ls->alpha);CHKERRQ(info);
  info=TaoPrintDouble(tao," beta=%g ",ls->beta);CHKERRQ(info);
  info=TaoPrintDouble(tao,"sigma=%g ",ls->sigma);CHKERRQ(info);
  info=TaoPrintDouble(tao,"minstep=%g,",ls->minimumStep);CHKERRQ(info);
  info=TaoPrintInt(tao,"memsize=%d\n",ls->memorySize);CHKERRQ(info);
  
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_PreArmijo"
static int TaoApply_PreArmijo(TAO_SOLVER tao, TAO_ARMIJO *ls, 
			      double f, double step, 
                              double *ref, TaoInt *idx, TaoInt *info2) 
{
  int info;
  TaoInt i;

  TaoFunctionBegin;

  *info2 = 0;

  // Check linesearch parameters
  if (step < 0) {
    info = PetscInfo1(tao, "TaoApply_Armijo:Line search error: step (%g) < 0\n", step); CHKERRQ(info);
    *info2 = -1; 
  } 

  if (ls->alpha < 1) {
    info = PetscInfo1(tao,"TaoApply_Armijo:Line search error: alpha (%g) < 1\n", ls->alpha); CHKERRQ(info);
    *info2 = -2; 
  } 
  
  if ((ls->beta <= 0) || (ls->beta >= 1)) {
    info = PetscInfo1(tao,"TaoApply_Armijo:Line search error: beta (%g) invalid\n", ls->beta); CHKERRQ(info);
    *info2 = -3; 
  } 
  
  if ((ls->beta_inf <= 0) || (ls->beta_inf >= 1)) {
    info = PetscInfo1(tao,"TaoApply_Armijo:Line search error: beta_inf (%g) invalid\n", ls->beta_inf); CHKERRQ(info);
    *info2 = -4; 
  } 

  if ((ls->sigma <= 0) || (ls->sigma >= 0.5)) {
    info = PetscInfo1(tao,"TaoApply_Armijo:Line search error: sigma (%g) invalid\n", ls->sigma); CHKERRQ(info);
    *info2 = -5; 
  } 
  
  if (ls->minimumStep <= 0) {
    info = PetscInfo1(tao,"TaoApply_Armijo:Line search error: minimum_step (%g) <= 0\n", ls->minimumStep); CHKERRQ(info);
    *info2 = -6; 
  } 
  
  if (ls->memorySize < 1) {
    info = PetscInfo1(tao,"TaoApply_Armijo:Line search error: memory_size (%d) < 1\n", ls->memorySize); CHKERRQ(info);
    *info2 = -7; 
  } 
  
  if ((ls->referencePolicy != REFERENCE_MAX) &&
      (ls->referencePolicy != REFERENCE_AVE) &&
      (ls->referencePolicy != REFERENCE_MEAN)) {
    info = PetscInfo(tao,"TaoApply_Armijo:Line search error: reference_policy invalid\n"); CHKERRQ(info);
    *info2 = -8; 
  } 
  
  if ((ls->replacementPolicy != REPLACE_FIFO) && 
      (ls->replacementPolicy != REPLACE_MRU)) {
    info = PetscInfo(tao,"TaoApply_Armijo:Line search error: replacement_policy invalid\n"); CHKERRQ(info);
    *info2 = -9; 
  }
  
  if (TaoInfOrNaN(f)) {
    info = PetscInfo(tao,"TaoApply_Armijo:Line search error: initial function inf or nan\n"); CHKERRQ(info);
    *info2 = -10; 
  }

  if (*info2) {
    TaoFunctionReturn(0);
  }

  // Check to see of the memory has been allocated.  If not, allocate
  // the historical array and populate it with the initial function
  // values.
  if (ls->memory == TAO_NULL) {
    info = TaoMalloc(sizeof(double)*ls->memorySize, &ls->memory ); CHKERRQ(info);
    
    info = PetscLogObjectMemory(tao, sizeof(double)*ls->memorySize); CHKERRQ(info);
  }

  if (tao->iter == 0) {
    for (i = 0; i < ls->memorySize; i++) {
      ls->memory[i] = ls->alpha*f;
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
  } 
  else if (ls->referencePolicy == REFERENCE_MEAN) {
    *ref = TaoMin(*ref, 0.5*(ls->lastReference + ls->memory[ls->current]));
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_PostArmijo"
static int TaoApply_PostArmijo(TAO_SOLVER tao, TAO_ARMIJO *ls, 
                               double f, double step,
			       double ref, TaoInt idx, TaoInt *info2) 
{
  int info;
  TaoFunctionBegin;

  *info2 = 0;

  // Check termination
  if (step < ls->minimumStep) {
    info = PetscInfo(tao, "TaoApply_Armijo:Step is at lower bound.\n"); CHKERRQ(info);
    *info2 = 1;
  }

  if (TaoInfOrNaN(f)) {
    info = PetscInfo(tao, "TaoApply_Armijo:Function is inf or nan.\n"); CHKERRQ(info);
    *info2 = 2;
  }

  if (*info2) {
    TaoFunctionReturn(0);
  }

  // Successful termination, update memory
  ls->lastReference = ref;
  if (ls->replacementPolicy == REPLACE_FIFO) {
    ls->memory[ls->current++] = f;
    if (ls->current >= ls->memorySize) {
      ls->current = 0;
    }
  } 
  else {
    ls->current = idx;
    ls->memory[idx] = f;
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_Armijo"
/* @ TaoApply_Armijo - This routine performs a linesearch. It
   backtracks until the (nonmonotone) Armijo conditions are satisfied.

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

static int TaoApply_Armijo(TAO_SOLVER tao, TaoVec *X, TaoVec *G, TaoVec *S,
                           TaoVec *W, double *f, double *f_full, double *step,
			   TaoInt *info2, void *ctx)
{
  TAO_ARMIJO *ls = (TAO_ARMIJO *)ctx;

  const double beta = ls->beta;
  const double beta_inf = ls->beta_inf;

  double fact, ref, t, gdx;
  TaoInt idx;
  int info;

  TaoFunctionBegin;
  info = TaoApply_PreArmijo(tao, ls, *f, *step, &ref, &idx, info2);

#if defined(PETSC_USE_COMPLEX)
  info = G->Dot(S,&cgdx);CHKERRQ(info); gdx = TaoReal(cgdx);
#else
  info = G->Dot(S,&gdx);CHKERRQ(info);
#endif

  if (TaoInfOrNaN(gdx)) {
    info = PetscInfo(tao,"TaoApply_Armijo:Line search error: gdx is inf or nan\n"); CHKERRQ(info);
    *info2 = -11;
  }

  if (gdx >= 0.0) {
    info = PetscInfo(tao,"TaoApply_LineSearch:Search direction not a descent direction\n"); CHKERRQ(info);
    *info2 = 12;
  }
  
  if (*info2) {
    TaoFunctionReturn(0);
  }

  fact = ls->sigma * gdx;
  t = *step;
  tao->new_search=TAO_TRUE;
  while (t >= ls->minimumStep) {
    // Calculate iterate
    info = W->Waxpby(1.0, X, t, S); CHKERRQ(info);

    // Calculate function at new iterate
    tao->current_step=t;
    info = TaoComputeMeritFunction(tao, W, f); CHKERRQ(info);
    tao->new_search=TAO_FALSE;
    if (*step == t) {
      *f_full = *f;
    }

    if (TaoInfOrNaN(*f)) {
      t *= beta_inf;
    }
    else {
      // Check descent condition
      if (*f <= ref + t*fact) {
        break;
      }

      t *= beta;
    }
  }

  info = TaoApply_PostArmijo(tao, ls, *f, t, ref, idx, info2);

  // Update iterate and compute gradient
  *step = t;
  info = X->CopyFrom(W); CHKERRQ(info);
  tao->current_step=t;
  info = TaoComputeMeritGradient(tao, X, G); CHKERRQ(info);

  // Finish computations
  info = PetscInfo1(tao, "TaoApply_Armijo:step = %10.4f\n",*step); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_NDArmijo"
/* @ TaoApply_NDArmijo - This routine performs a linesearch. It
   backtracks until the (nonmonotone) Armijo conditions are satisfied.
   This is modified for a nondifferentiable function.

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

static int TaoApply_NDArmijo(TAO_SOLVER tao, TaoVec *X, TaoVec *G, TaoVec *S,
                             TaoVec *W, double *f, double *f_full, double *step,
			     TaoInt *info2, void *ctx)
{
  TAO_ARMIJO *ls = (TAO_ARMIJO *)ctx;

  const double fact = ls->sigma;
  const double beta = ls->beta;
  const double beta_inf = ls->beta_inf;

  double ref, t;
  TaoInt idx;
  int info;

  TaoFunctionBegin;

  info = TaoApply_PreArmijo(tao, ls, *f, *step, &ref, &idx, info2);
  if (*info2) {
    TaoFunctionReturn(0);
  }

  t = *step;
  tao->new_search=TAO_TRUE;
  while (t >= ls->minimumStep) {
    // Calculate iterate
    info = W->Waxpby(1.0, X, t, S); CHKERRQ(info);

    // Calculate function at new iterate
    tao->current_step=t;
    info = TaoComputeMeritFunction(tao, W, f); CHKERRQ(info);
    tao->new_search=TAO_FALSE;
    if (*step == t) {
      *f_full = *f;
    }

    if (TaoInfOrNaN(*f)) {
      t *= beta_inf;
    }
    else { 
      // Check descent condition
      if (*f <= (1 - fact*t)*ref) {
        break;
      }

      t *= beta;
    }
  }

  info = TaoApply_PostArmijo(tao, ls, *f, t, ref, idx, info2);

  // Update iterate and compute gradient
  *step = t;
  info = X->CopyFrom(W); CHKERRQ(info);
  tao->current_step=t;
  info = TaoComputeMeritGradient(tao, X, G); CHKERRQ(info);

  // Finish computations
  info = PetscInfo1(tao, "TaoApply_NDArmijo:step = %10.4f\n",*step); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCreateArmijoLineSearch"
/*@C
   TaoCreateArmijoLineSearch - Create a non-monotone linesearch

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
   unconstrained optimization.  This routine is used within the following
   TAO solvers: infeasible semismooth with linesearch (tao_ssils).

   Level: developer

.keywords: TAO_SOLVER, linesearch
@*/
int TaoCreateArmijoLineSearch(TAO_SOLVER tao)
{
  TAO_ARMIJO *ls;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_ARMIJO, &ls);CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_ARMIJO)); CHKERRQ(info);

  ls->memory = TAO_NULL;
  ls->alpha = 1.0;
  ls->beta = 0.5;
  ls->beta_inf = 0.5;
  ls->sigma = 1e-4;
  ls->minimumStep = TAO_EPSILON;
  ls->memorySize = 1;
  ls->referencePolicy = REFERENCE_MAX;
  ls->replacementPolicy = REPLACE_MRU;

  info = TaoSetLineSearch(tao,0,
			  TaoSetOptions_Armijo,
			  TaoApply_Armijo,
			  TaoView_Armijo,
			  TaoDestroy_Armijo,
			  (void *) ls);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCreateNDArmijoLineSearch"
/*@C
   TaoCreateNDArmijoLineSearch - Create a non-monotone linesearch for a 
     nondifferentiable function

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
   unconstrained optimization.  This routine is used within the following
   TAO solvers: infeasible semismooth with linesearch (tao_ssils).

   Level: developer

.keywords: TAO_SOLVER, linesearch
@*/
int TaoCreateNDArmijoLineSearch(TAO_SOLVER tao)
{
  TAO_ARMIJO *ls;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_ARMIJO, &ls);CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_ARMIJO)); CHKERRQ(info);

  ls->memory = TAO_NULL;
  ls->alpha = 1.0;
  ls->beta = 0.5;
  ls->beta_inf = 0.5;
  ls->sigma = 1e-4;
  ls->minimumStep = TAO_EPSILON;
  ls->memorySize = 1;
  ls->referencePolicy = REFERENCE_MAX;
  ls->replacementPolicy = REPLACE_MRU;

  info = TaoSetLineSearch(tao,0,
			  TaoSetOptions_Armijo,
			  TaoApply_NDArmijo,
			  TaoView_Armijo,
			  TaoDestroy_Armijo,
			  (void *) ls);CHKERRQ(info);

  TaoFunctionReturn(0);
}

