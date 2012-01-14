#include "src/complementarity/impls/ssls/ssls.h"

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_SSFLS"
static int TaoSolve_SSFLS(TAO_SOLVER tao, void *solver)
{
  TAO_SSLS *ssls = (TAO_SSLS *)solver;
  //  TaoLinearSolver *lsolver;
  TaoVec *x, *l, *u, *ff, *dpsi, *d, *w;
  TaoMat *J;
  double psi, psi_full, ndpsi, normd, innerd, t=0;
  double delta, rho;
  int iter=0, info;
  TaoTerminateReason reason;
  TaoTruth flag;

  TaoFunctionBegin;

  // Assume that Setup has been called!
  // Set the structure for the Jacobian and create a linear solver.
 
  delta = ssls->delta;
  rho = ssls->rho;

  info = TaoGetSolution(tao, &x); CHKERRQ(info);
  l=ssls->xl;
  u=ssls->xu;
  info = TaoEvaluateVariableBounds(tao,l,u); CHKERRQ(info);
  info = x->Median(l,x,u); CHKERRQ(info);
  info = TaoGetJacobian(tao, &J); CHKERRQ(info);

  ff = ssls->ff;
  dpsi = ssls->dpsi;
  d = ssls->d;
  w = ssls->w;

  info = x->PointwiseMaximum(x, l); CHKERRQ(info);
  info = x->PointwiseMinimum(x, u); CHKERRQ(info);
  info = TaoSetMeritFunction(tao, Tao_SSLS_Function, Tao_SSLS_FunctionGradient,
			     TAO_NULL, TAO_NULL, TAO_NULL, ssls); CHKERRQ(info);

  // Calculate the function value and fischer function value at the 
  // current iterate
  info = TaoComputeMeritFunctionGradient(tao, x, &psi, dpsi); CHKERRQ(info);
  info = dpsi->Norm2(&ndpsi);

  while (1) {
    info=PetscInfo3(tao, "TaoSolve_SSFLS: %d: merit: %5.4e, ndpsi: %5.4e\n",
		       iter, ssls->merit, ndpsi);CHKERRQ(info);

    // Check the termination criteria
    info = TaoMonitor(tao,iter++,ssls->merit,ndpsi,0.0,t,&reason); 
           CHKERRQ(info);
    if (reason!=TAO_CONTINUE_ITERATING) break;

    // Calculate direction.  (Really negative of newton direction.  Therefore,
    // rest of the code uses -d.)
    info = TaoPreLinearSolve(tao, J); CHKERRQ(info);
    info = TaoLinearSolve(tao, J, ff, d, &flag); CHKERRQ(info);
    
    info = w->CopyFrom(d); CHKERRQ(info);
    info = w->Negate(); CHKERRQ(info);
    info = w->BoundGradientProjection(w,l, x, u);

    info = w->Norm2(&normd); CHKERRQ(info);
    info = w->Dot(dpsi, &innerd); CHKERRQ(info);

    // Make sure that we have a descent direction
    if (innerd >= -delta*pow(normd, rho)) {
      info = PetscInfo1(tao, "TaoSolve_SSFLS: %d: newton direction not descent\n", iter); CHKERRQ(info);
      info = d->CopyFrom(dpsi); CHKERRQ(info);
      info = w->Dot(dpsi, &innerd); CHKERRQ(info);
    }
    info = d->Negate(); CHKERRQ(info);
    innerd = -innerd;

    t = 1;
    info = TaoLineSearchApply(tao, x, dpsi, d, w, 
                              &psi, &psi_full, &t, &tao->lsflag); CHKERRQ(info);
    info = dpsi->Norm2(&ndpsi);
  }
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_SSFLS"
int TaoCreate_SSFLS(TAO_SOLVER tao)
{
  TAO_SSLS *ssls;
  int        info;

  TaoFunctionBegin;

  info = TaoNew(TAO_SSLS,&ssls); CHKERRQ(info);
  info = PetscLogObjectMemory(tao, sizeof(TAO_SSLS)); CHKERRQ(info);

  ssls->delta = 1e-10;
  ssls->rho = 2.1;

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_SSFLS,(void*)ssls); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_SSLS,TaoSetDown_SSLS); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetOptions_SSLS); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_SSLS); CHKERRQ(info);

  info = TaoCreateProjectedArmijoLineSearch(tao); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,2000); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,4000); CHKERRQ(info);

  info = TaoSetTolerances(tao,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,1.0e-16,0.0,0.0); CHKERRQ(info);
  info = TaoSetFunctionLowerBound(tao,1.0e-8); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END

