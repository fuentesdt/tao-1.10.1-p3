/*$Id$*/

#include "cg.h"

#define CG_FletcherReeves       0
#define CG_PolakRibiere         1
#define CG_PolakRibierePlus     2
#define CG_HestenesStiefel      3
#define CG_DaiYuan              4
#define CG_Types                5

static const char *CG_Table[64] = {
  "fr", "pr", "prp", "hs", "dy"
};

#define TAO_ZER_SAFEGUARD	1e-8
#define TAO_INF_SAFEGUARD	1e+8

#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_CG"
static int TaoSolve_CG(TAO_SOLVER tao, void *solver)
{
  TAO_CG *cg = (TAO_CG *)solver;
  TaoVec *X, *Xm1 = cg->X2;
  TaoVec *G = cg->G1, *Gm1 = cg->G2;
  TaoVec *D = cg->D, *W = cg->W;

  TaoTerminateReason reason;

  double f, f_full, fm1, gnorm, gnorm2, gnorm2m1, ginner, gd, gm1d;
  double beta, delta, step = 1.0;

  int info = 0;
  TaoInt status = 0, iter=0;

  TaoFunctionBegin;

  // Get vectors we will need
  info = TaoGetSolution(tao, &X); CHKERRQ(info);

  // Check convergence criteria
  info = TaoComputeFunctionGradient(tao, X, &f, G); CHKERRQ(info);
  info = G->Norm2(&gnorm); CHKERRQ(info);
  if (TaoInfOrNaN(f) || TaoInfOrNaN(gnorm)) {
    SETERRQ(1, "User provided compute function generated Inf or NaN");
  }

  info = TaoMonitor(tao, iter, f, gnorm, 0.0, step, &reason); CHKERRQ(info);
  if (reason != TAO_CONTINUE_ITERATING) {
    TaoFunctionReturn(0);
  }

  // Have not converged; initialize variables
  info = D->ScaleCopyFrom(-1.0, G); CHKERRQ(info);
  gnorm2 = gnorm*gnorm;

  // Set initial scaling for the function
  if (f != 0.0) {
    delta = 2.0 * TaoAbsDouble(f) / gnorm2;
    delta = TaoMax(delta, cg->delta_min);
    delta = TaoMin(delta, cg->delta_max);
  }
  else {
    delta = 2.0 / gnorm2;
    delta = TaoMax(delta, cg->delta_min);
    delta = TaoMin(delta, cg->delta_max);
  }

  // Set counter for gradient/reset steps
  cg->grad = 0;
  cg->reset = 0;

  while (1) {
    // Save the current gradient information
    fm1 = f;
    gnorm2m1 = gnorm2;
    info = Xm1->CopyFrom(X); CHKERRQ(info);
    info = Gm1->CopyFrom(G); CHKERRQ(info);

    info = D->Dot(G, &gd); CHKERRQ(info);
    if ((gd >= 0) || TaoInfOrNaN(gd)) {
      // Step is not descent or direction generated not a number
      // Use steepest descent direction
      ++cg->grad;

      if (f != 0.0) {
        delta = 2.0 * TaoAbsDouble(f) / gnorm2;
        delta = TaoMax(delta, cg->delta_min);
        delta = TaoMin(delta, cg->delta_max);
      }
      else {
        delta = 2.0 / gnorm2;
        delta = TaoMax(delta, cg->delta_min);
        delta = TaoMin(delta, cg->delta_max);
      }

      info = D->ScaleCopyFrom(-1.0, G); CHKERRQ(info);

      // Gradient step cannot include not a number; this test is not needed.
      // info = D->Norm2(&dnorm); CHKERRQ(info);
      // if (TaoInfOrNaN(dnorm)) {
      //   SETERRQ(1, "Direction generated Not-a-Number");
      // }
    }

    // Search direction for improving point
    step = delta;
    info = TaoLineSearchApply(tao, X, G, D, W, &f, &f_full, &step, &status); CHKERRQ(info);

    if (status) {
      // Linesearch failed
      // Reset factors and use scaled gradient step
      ++cg->reset;

      f = fm1;
      gnorm2 = gnorm2m1;
      info = X->CopyFrom(Xm1); CHKERRQ(info);
      info = G->CopyFrom(Gm1); CHKERRQ(info);

      if (f != 0.0) {
        delta = 2.0 * TaoAbsDouble(f) / gnorm2;
        delta = TaoMax(delta, cg->delta_min);
        delta = TaoMin(delta, cg->delta_max);
      }
      else {
        delta = 2.0 / gnorm2;
        delta = TaoMax(delta, cg->delta_min);
        delta = TaoMin(delta, cg->delta_max);
      }

      info = D->ScaleCopyFrom(-1.0, G); CHKERRQ(info);

      // Gradient step cannot include not a number; this test is not needed.
      // info = D->Norm2(&dnorm); CHKERRQ(info);
      // if (TaoInfOrNaN(dnorm)) {
      //   SETERRQ(1, "Direction generated Not-a-Number");
      // }

      // This may be incorrect; linesearch has values for stepmax and stepmin
      // that should be reset.
      step = delta;
      info = TaoLineSearchApply(tao, X, G, D, W, &f, &f_full, &step, &status); CHKERRQ(info);

      if (status) {
        // Linesearch failed,
        // Switch to unscaled gradient

        f = fm1;
        gnorm2 = gnorm2m1;
        info = X->CopyFrom(Xm1); CHKERRQ(info);
        info = G->CopyFrom(Gm1); CHKERRQ(info);

        delta = 1.0;

        info = D->ScaleCopyFrom(-1.0, G); CHKERRQ(info);

        // Gradient step cannot include not a number; this test is not needed.
        // info = D->Norm2(&dnorm); CHKERRQ(info);
        // if (TaoInfOrNaN(dnorm)) {
        //   SETERRQ(1, "Direction generated Not-a-Number");
        // }

        // This may be incorrect; linesearch has values for stepmax and stepmin
        // that should be reset.
        step = delta;
        info = TaoLineSearchApply(tao, X, G, D, W, &f, &f_full, &step, &status); CHKERRQ(info);
        if (status) {
          // Steepest descent direction did not produce a new value
          // Stop here

          f = fm1;
          gnorm2 = gnorm2m1;
          info = X->CopyFrom(Xm1); CHKERRQ(info);
          info = G->CopyFrom(Gm1); CHKERRQ(info);
          step = 0.0;
        }
      }
    }

    // Check for termination
    info = G->Norm2(&gnorm); CHKERRQ(info);
    if (TaoInfOrNaN(f) || TaoInfOrNaN(gnorm)) {
      SETERRQ(1, "User provided compute function generated Inf or NaN");
    }
    gnorm2 = gnorm*gnorm;

    // Check for termination
    info = TaoMonitor(tao, ++iter, f, gnorm, 0.0, step, &reason); CHKERRQ(info);
    if (reason != TAO_CONTINUE_ITERATING) {
      break;
    }

    // Check for restart condition
    info = G->Dot(Gm1, &ginner); CHKERRQ(info);
    if (fabs(ginner) >= cg->eta * gnorm2) {
      // Gradients far from orthogonal; use steepest descent direction
      beta = 0.0;
    }
    else {
      // Gradients close to orthogonal; use conjugate gradient formula

      switch(cg->cg_type) {
      case CG_FletcherReeves:
	beta = gnorm2 / gnorm2m1;
	break;

      case CG_PolakRibiere:
	beta = (gnorm2 - ginner) / gnorm2m1;
	break;

      case CG_PolakRibierePlus:
	beta = TaoMax((gnorm2 - ginner) / gnorm2m1, 0.0);
        break;

      case CG_HestenesStiefel:
	info = G->Dot(D, &gd); CHKERRQ(info);
	info = Gm1->Dot(D, &gm1d); CHKERRQ(info);
	beta = (gnorm2 - ginner) / (gd - gm1d); 
	break;

      case CG_DaiYuan:
	info = G->Dot(D, &gd); CHKERRQ(info);
	info = Gm1->Dot(D, &gm1d); CHKERRQ(info);
	beta = gnorm2 / (gd - gm1d); 
	break;

      default:
        beta = 0.0;
        break;
      }
    }

    // Compute the direction
    info = D->Axpby(-1.0, G, beta); CHKERRQ(info);

    // Update initial steplength choice
    delta = 1.0;
    delta = TaoMax(delta, cg->delta_min);
    delta = TaoMin(delta, cg->delta_max);
  }
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_CG"
static int TaoSetUp_CG(TAO_SOLVER tao, void *solver)
{
  TAO_CG *cg = (TAO_CG *)solver;
  TaoVec *X;
  int info;

  TaoFunctionBegin;
  
  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = X->Clone(&cg->X2); CHKERRQ(info);
  info = X->Clone(&cg->G1); CHKERRQ(info);
  info = X->Clone(&cg->G2); CHKERRQ(info);
  info = X->Clone(&cg->D); CHKERRQ(info);
  info = X->Clone(&cg->W); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao, cg->G1); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, cg->D); CHKERRQ(info);

  info = TaoCheckFG(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_CG"
static int TaoDestroy_CG(TAO_SOLVER tao, void *solver)
{
  TAO_CG *cg = (TAO_CG *)solver;
  int info;

  TaoFunctionBegin;

  info = TaoVecDestroy(cg->X2); CHKERRQ(info);
  info = TaoVecDestroy(cg->G1); CHKERRQ(info);
  info = TaoVecDestroy(cg->G2); CHKERRQ(info);
  info = TaoVecDestroy(cg->D); CHKERRQ(info);
  info = TaoVecDestroy(cg->W); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao, 0); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, 0); CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_CG"
static int TaoSetOptions_CG(TAO_SOLVER tao, void *solver)
{
  TAO_CG *cg = (TAO_CG *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoOptionsHead("Nonlinear Conjugate Gradient method for unconstrained optimization"); CHKERRQ(info);

  info = TaoOptionDouble("-tao_cg_eta", "restart tolerance", "", cg->eta, &cg->eta, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_cg_type", "cg formula", "", CG_Table, CG_Types, CG_Table[cg->cg_type], &cg->cg_type, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_cg_delta_min", "minimum delta value", "", cg->delta_min, &cg->delta_min, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_cg_delta_max", "maximum delta value", "", cg->delta_max, &cg->delta_max, 0); CHKERRQ(info);

  info = TaoLineSearchSetFromOptions(tao); CHKERRQ(info);
  info = TaoOptionsTail(); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_CG"
static int TaoView_CG(TAO_SOLVER tao, void *solver)
{
  TAO_CG   *cg = (TAO_CG *)solver;
  int      info;

  TaoFunctionBegin;
  info = TaoPrintInt(tao, "  CG Type: %d\n", cg->cg_type); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Gradient steps: %d\n", cg->grad); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Reset steps: %d\n", cg->reset); CHKERRQ(info);
  info = TaoLineSearchView(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_CG"
int TaoCreate_CG(TAO_SOLVER tao)
{
  TAO_CG *cg;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_CG, &cg); CHKERRQ(info);
  info = PetscLogObjectMemory(tao, sizeof(TAO_CG)); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao, TaoSolve_CG, (void *)cg); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao, TaoSetUp_CG, TaoDestroy_CG); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao, TaoSetOptions_CG); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao, TaoView_CG); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao, 2000); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao, 4000); CHKERRQ(info);
  info = TaoSetTolerances(tao, 1e-4, 1e-4, 0, 0); CHKERRQ(info);

  cg->eta = 0.1;
  cg->delta_min = 1e-7;
  cg->delta_max = 100;

  cg->cg_type = CG_PolakRibierePlus;

  // Note: nondefault values should be used for nonlinear conjugate gradient 
  // method.  In particular, gtol should be less that 0.5; the value used in 
  // Nocedal and Wright is 0.10.  We use the default values for the 
  // linesearch because it seems to work better.
  info = TaoCreateMoreThuenteLineSearch(tao, 0, 0); CHKERRQ(info);
  TaoFunctionReturn(0);
}
EXTERN_C_END
