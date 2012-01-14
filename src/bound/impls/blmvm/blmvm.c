/*$Id$*/

#include "blmvm.h"

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_BLMVM"
static int TaoSolve_BLMVM(TAO_SOLVER tao, void *solver)
{
  TAO_BLMVM *blm = (TAO_BLMVM *)solver;
  TaoVec *X, *XL = blm->XL, *XU = blm->XU;
  TaoVec *G = blm->G, *GP = blm->GP, *D = blm->D;
  TaoVec *Xold, *Gold;
  TaoLMVMMat *M = blm->M;

  TaoTerminateReason reason;
  TaoTruth success;

  double f, f_full, fold, gdx, gnorm;
  double step = 1.0;
  int info;
  TaoInt iter = 0, status = 0;
  
  TaoFunctionBegin;
  
  // Get vectors we will need
  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = X->Clone(&Xold); CHKERRQ(info);
  info = X->Clone(&Gold); CHKERRQ(info);

  info = TaoEvaluateVariableBounds(tao, XL, XU); CHKERRQ(info);
  
  // Project initial point onto bounds
  info = X->Median(XL, X, XU); CHKERRQ(info);

  // Check convergence criteria
  info = TaoComputeMeritFunctionGradient(tao, X, &f, G); CHKERRQ(info);
  info = GP->BoundGradientProjection(G, XL, X, XU); CHKERRQ(info);
  info = GP->Norm2(&gnorm); CHKERRQ(info);
  if ((f != f) || (gnorm != gnorm)) {
    SETERRQ(1, "User provided compute function generated Not-a-Number");
  }

  info = TaoMonitor(tao, iter, f, gnorm, 0.0, step, &reason); CHKERRQ(info);
  if (reason != TAO_CONTINUE_ITERATING) {
    TaoFunctionReturn(0);
  }

  // Set initial scaling for the function
  if (f != 0.0) {
    info = M->SetDelta(2.0 * TaoAbsDouble(f) / (gnorm*gnorm)); CHKERRQ(info);
  }
  else {
    info = M->SetDelta(2.0 / (gnorm*gnorm)); CHKERRQ(info);
  }
 
  // Set counter for gradient/reset steps
  blm->grad = 0;
  blm->reset = 0;

  // Have not converged; continue with Newton method
  while (reason == TAO_CONTINUE_ITERATING) {
    
    // Compute direction
    info = M->Update(X, GP); CHKERRQ(info);
    info = M->Solve(G, D, &success); CHKERRQ(info);
    info = GP->BoundGradientProjection(D, XL, X, XU); CHKERRQ(info);

    // Check for success (descent direction)
    info = GP->Dot(G, &gdx); CHKERRQ(info);
    if (gdx <= 0) {
      // Step is not descent or solve was not successful
      // Use steepest descent direction (scaled)
      ++blm->grad;

      if (f != 0.0) {
        info = M->SetDelta(2.0 * TaoAbsDouble(f) / (gnorm*gnorm)); CHKERRQ(info);
      }
      else {
        info = M->SetDelta(2.0 / (gnorm*gnorm)); CHKERRQ(info);
      }
      info = M->Reset(); CHKERRQ(info);
      info = M->Update(X, G); CHKERRQ(info);
      info = M->Solve(G, D, &success); CHKERRQ(info);
    } 
    info = D->Negate(); CHKERRQ(info);

    // Perform the linesearch
    fold = f;
    info = Xold->CopyFrom(X); CHKERRQ(info);
    info = Gold->CopyFrom(G); CHKERRQ(info);

    step = 1.0;
    info = TaoLineSearchApply(tao, X, G, D, GP, &f, &f_full, &step, &status); CHKERRQ(info);

    if (status) {
      // Linesearch failed
      // Reset factors and use scaled (projected) gradient step
      ++blm->reset;

      f = fold;
      info = X->CopyFrom(Xold); CHKERRQ(info);
      info = G->CopyFrom(Gold); CHKERRQ(info);

      if (f != 0.0) {
        info = M->SetDelta(2.0 * TaoAbsDouble(f) / (gnorm*gnorm)); CHKERRQ(info);
      }
      else {
        info = M->SetDelta(2.0 / (gnorm*gnorm)); CHKERRQ(info);
      }
      info = M->Reset(); CHKERRQ(info);
      info = M->Update(X, G); CHKERRQ(info);
      info = M->Solve(G, D, &success); CHKERRQ(info);
      info = D->Negate(); CHKERRQ(info);

      // This may be incorrect; linesearch has values fo stepmax and stepmin
      // that should be reset.
      step = 1.0;
      info = TaoLineSearchApply(tao, X, G, D, GP, &f, &f_full, &step, &status); CHKERRQ(info);

      if (status) {
        // Linesearch failed
        // Probably stop here
      }
    }

    // Check for termination
    info = GP->BoundGradientProjection(G, XL, X, XU); CHKERRQ(info);
    info = GP->Norm2(&gnorm); CHKERRQ(info);
    if ((f != f) || (gnorm != gnorm)) {
      SETERRQ(1, "User provided compute function generated Not-a-Number");
    }
    info = TaoMonitor(tao, ++iter, f, gnorm, 0.0, step, &reason); CHKERRQ(info);
  }

  info = TaoVecDestroy(Xold); CHKERRQ(info);
  info = TaoVecDestroy(Gold); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_BLMVM"
static int TaoSetUp_BLMVM(TAO_SOLVER tao, void *solver)
{
  TAO_BLMVM *blm = (TAO_BLMVM *)solver;
  TaoVec *X;
  int info;

  TaoFunctionBegin;

  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = X->Clone(&blm->XL); CHKERRQ(info);
  info = X->Clone(&blm->XU); CHKERRQ(info);
  info = X->Clone(&blm->D); CHKERRQ(info);
  info = X->Clone(&blm->G); CHKERRQ(info);
  info = X->Clone(&blm->GP); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao, blm->GP); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, blm->D); CHKERRQ(info);
  info = TaoSetVariableBounds(tao, blm->XL, blm->XU); CHKERRQ(info);
  
  // Create matrix for the limited memory approximation
  blm->M = new TaoLMVMMat(X);

  info = TaoCheckFG(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_BLMVM"
static int TaoSetDown_BLMVM(TAO_SOLVER tao, void *solver)
{
  TAO_BLMVM *blm = (TAO_BLMVM *)solver;
  int info;

  TaoFunctionBegin;

  info=TaoVecDestroy(blm->XL); CHKERRQ(info);
  info=TaoVecDestroy(blm->XU); CHKERRQ(info);
  info=TaoVecDestroy(blm->G); CHKERRQ(info);
  info=TaoVecDestroy(blm->GP); CHKERRQ(info);
  info=TaoVecDestroy(blm->D); CHKERRQ(info);
  info=TaoMatDestroy(blm->M); CHKERRQ(info);
  
  info = TaoSetLagrangianGradientVector(tao, 0); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, 0); CHKERRQ(info);
  info = TaoSetVariableBounds(tao, 0, 0); CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_BLMVM"
static int TaoSetOptions_BLMVM(TAO_SOLVER tao, void *solver)
{
  int info;

  TaoFunctionBegin;
  info = TaoOptionsHead("Limited-memory variable-metric method for bound constrained optimization"); CHKERRQ(info);
  info = TaoLineSearchSetFromOptions(tao);CHKERRQ(info);
  info = TaoOptionsTail();CHKERRQ(info);
  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_BLMVM"
static int TaoView_BLMVM(TAO_SOLVER tao, void *solver)
{
  TAO_BLMVM *blm = (TAO_BLMVM *) solver;
  int info;

  TaoFunctionBegin;
  info = TaoPrintInt(tao, "  Rejected matrix updates: %d\n", blm->M->GetRejects()); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Gradient steps: %d\n", blm->grad); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Reset steps: %d\n", blm->reset); CHKERRQ(info);
  info = TaoLineSearchView(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetDualVariables_BLMVM" 
static int TaoGetDualVariables_BLMVM(TAO_SOLVER tao, TaoVec *DXL, TaoVec* DXU, void *solver)
{
  TAO_BLMVM *blm = (TAO_BLMVM *) solver;
  TaoVec *G = blm->G, *GP = blm->GP;
  int info;

  TaoFunctionBegin;
  info = DXL->Waxpby(-1.0, G, 1.0, GP); CHKERRQ(info);
  info = DXU->SetToZero(); CHKERRQ(info);
  info = DXL->PointwiseMaximum(DXL, DXU); CHKERRQ(info);

  info = DXU->Waxpby(-1.0, GP, 1.0, G); CHKERRQ(info);
  info = DXU->Axpy(1.0, DXL); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_BLMVM"
int TaoCreate_BLMVM(TAO_SOLVER tao)
{
  TAO_BLMVM *blm;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_BLMVM, &blm); CHKERRQ(info);
  info = PetscLogObjectMemory(tao, sizeof(TAO_BLMVM)); CHKERRQ(info);

  info = TaoSetTaoSolveRoutine(tao, TaoSolve_BLMVM, (void*)blm); CHKERRQ(info);
  info = TaoSetTaoSetUpDownRoutines(tao, TaoSetUp_BLMVM, TaoSetDown_BLMVM); CHKERRQ(info);
  info = TaoSetTaoOptionsRoutine(tao, TaoSetOptions_BLMVM); CHKERRQ(info);
  info = TaoSetTaoViewRoutine(tao, TaoView_BLMVM); CHKERRQ(info);
  info = TaoSetTaoDualVariablesRoutine(tao, TaoGetDualVariables_BLMVM); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao, 2000); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao, 4000); CHKERRQ(info);
  info = TaoSetTolerances(tao, 1e-4, 1e-4, 0, 0); CHKERRQ(info);

  info = TaoCreateMoreThuenteBoundLineSearch(tao, 0, 0); CHKERRQ(info);
  TaoFunctionReturn(0);
}
EXTERN_C_END

