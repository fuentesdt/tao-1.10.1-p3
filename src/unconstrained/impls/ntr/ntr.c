/*$Id$*/

#include "ntr.h"

#ifdef TAO_USE_PETSC
#include "petscksp.h"
#include "petscpc.h"
#include "src/petsctao/linearsolver/taolinearsolver_petsc.h"
#include "src/petsctao/vector/taovec_petsc.h"

#include "private/kspimpl.h"
#include "private/pcimpl.h"

#define NTR_KSP_NASH    0
#define NTR_KSP_STCG    1
#define NTR_KSP_GLTR    2
#define NTR_KSP_TYPES   3

#define NTR_PC_NONE	0
#define NTR_PC_AHESS    1
#define NTR_PC_BFGS     2
#define NTR_PC_PETSC    3
#define NTR_PC_TYPES    4

#define BFGS_SCALE_AHESS   0
#define BFGS_SCALE_BFGS    1
#define BFGS_SCALE_TYPES   2

#define NTR_INIT_CONSTANT	  0
#define NTR_INIT_DIRECTION	  1
#define NTR_INIT_INTERPOLATION	  2
#define NTR_INIT_TYPES		  3

#define NTR_UPDATE_REDUCTION      0
#define NTR_UPDATE_INTERPOLATION  1
#define NTR_UPDATE_TYPES          2

static const char *NTR_KSP[64] = {
  "nash", "stcg", "gltr"
};

static const char *NTR_PC[64] = {
  "none", "ahess", "bfgs", "petsc"
};

static const char *BFGS_SCALE[64] = {
  "ahess", "bfgs"
};

static const char *NTR_INIT[64] = {
  "constant", "direction", "interpolation"
};

static const char *NTR_UPDATE[64] = {
  "reduction", "interpolation"
};

// Routine for BFGS preconditioner

#undef __FUNCT__  
#define __FUNCT__ "bfgs_apply"
static PetscErrorCode bfgs_apply(PC pc, Vec xin, Vec xout)
{
  TaoLMVMMat *M;
  TaoVecPetsc Xin(xin);
  TaoVecPetsc Xout(xout);
  TaoTruth info2;
  int info;

  PetscFunctionBegin;
  info = PCShellGetContext(pc,(void**)&M); CHKERRQ(info);
  info = M->Solve(&Xin, &Xout, &info2); CHKERRQ(info);
  PetscFunctionReturn(0);
}

// TaoSolve_NTR - Implements Newton's Method with a trust region approach 
// for solving unconstrained minimization problems.  
//
// The basic algorithm is taken from MINPACK-2 (dstrn).
//
// TaoSolve_NTR computes a local minimizer of a twice differentiable function
// f by applying a trust region variant of Newton's method.  At each stage 
// of the algorithm, we use the prconditioned conjugate gradient method to
// determine an approximate minimizer of the quadratic equation
//
//      q(s) = <s, Hs + g>
//
// subject to the trust region constraint
//
//      || s ||_M <= radius,
//
// where radius is the trust region radius and M is a symmetric positive
// definite matrix (the preconditioner).  Here g is the gradient and H 
// is the Hessian matrix.
//
// Note:  TaoSolve_NTR MUST use the iterative solver KSPNASH, KSPSTCG,
//        or KSPGLTR.  Thus, we set KSPNASH, KSPSTCG, or KSPGLTR in this 
//        routine regardless of what the user may have previously specified.

#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_NTR"
static int TaoSolve_NTR(TAO_SOLVER tao, void *solver)
{
  TAO_NTR *tr = (TAO_NTR *)solver;
  TaoVec *X, *G = tr->G, *D = tr->D, *W = tr->W;
  TaoVec *Diag = tr->Diag;
  TaoMat *H;
  TaoLMVMMat *M = tr->M;

  TaoLinearSolver *tls;
  TaoLinearSolverPetsc *pls;

  KSP pksp;
  PC ppc;

  KSPConvergedReason ksp_reason;
  TaoTerminateReason reason;
  TaoTruth success;

  double fmin, ftrial, prered, actred, kappa, sigma, beta;
  double tau, tau_1, tau_2, tau_max, tau_min, max_radius;
  double f, gnorm;

  double delta;
  double radius, norm_d;
  int info;

  TaoInt iter = 0;
  TaoInt bfgsUpdates = 0;
  TaoInt needH;

  TaoInt i_max = 5;
  TaoInt j_max = 1;
  TaoInt i, j;

  TaoFunctionBegin;

  // Get the initial trust-region radius
  info = TaoGetInitialTrustRegionRadius(tao, &radius); CHKERRQ(info);
  if (radius < 0.0) {
    SETERRQ(1, "Initial radius negative");
  }

  // Modify the radius if it is too large or small
  radius = TaoMax(radius, tr->min_radius);
  radius = TaoMin(radius, tr->max_radius);

  // Get vectors we will need
  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = TaoGetHessian(tao, &H); CHKERRQ(info);

  if (NTR_PC_BFGS == tr->pc_type && !M) {
    tr->M = new TaoLMVMMat(X);
    M = tr->M;
  }

  // Check convergence criteria
  info = TaoComputeFunctionGradient(tao, X, &f, G); CHKERRQ(info);
  info = G->Norm2(&gnorm); CHKERRQ(info);
  if (TaoInfOrNaN(f) || TaoInfOrNaN(gnorm)) {
    SETERRQ(1, "User provided compute function generated Inf or NaN");
  }
  needH = 1;

  info = TaoMonitor(tao, iter, f, gnorm, 0.0, 1.0, &reason); CHKERRQ(info);
  if (reason != TAO_CONTINUE_ITERATING) {
    TaoFunctionReturn(0);
  }

  // Create vectors for the limited memory preconditioner
  if ((NTR_PC_BFGS == tr->pc_type) && 
      (BFGS_SCALE_BFGS != tr->bfgs_scale_type)) {
    if (!Diag) {
      info = X->Clone(&tr->Diag); CHKERRQ(info);
      Diag = tr->Diag;
    }
  }

  // Modify the linear solver to a conjugate gradient method
  info = TaoGetLinearSolver(tao, &tls); CHKERRQ(info);
  pls  = dynamic_cast <TaoLinearSolverPetsc *> (tls);

  pksp = pls->GetKSP();
  switch(tr->ksp_type) {
  case NTR_KSP_NASH:
    info = KSPSetType(pksp, KSPNASH); CHKERRQ(info);
    if (pksp->ops->setfromoptions) {
      (*pksp->ops->setfromoptions)(pksp);
    }
    break;

  case NTR_KSP_STCG:
    info = KSPSetType(pksp, KSPSTCG); CHKERRQ(info);
    if (pksp->ops->setfromoptions) {
      (*pksp->ops->setfromoptions)(pksp);
    }
    break;

  default:
    info = KSPSetType(pksp, KSPGLTR); CHKERRQ(info);
    if (pksp->ops->setfromoptions) {
      (*pksp->ops->setfromoptions)(pksp);
    }
    break;
  }

  // Modify the preconditioner to use the bfgs approximation
  info = KSPGetPC(pksp, &ppc); CHKERRQ(info);
  switch(tr->pc_type) {
  case NTR_PC_NONE:
    info = PCSetType(ppc, PCNONE); CHKERRQ(info);
    if (ppc->ops->setfromoptions) {
      (*ppc->ops->setfromoptions)(ppc);
    }
    break;
 
  case NTR_PC_AHESS:
    info = PCSetType(ppc, PCJACOBI); CHKERRQ(info);
    if (ppc->ops->setfromoptions) {
      (*ppc->ops->setfromoptions)(ppc);
    }
    info = PCJacobiSetUseAbs(ppc); CHKERRQ(info);
    break;

  case NTR_PC_BFGS:
    info = PCSetType(ppc, PCSHELL); CHKERRQ(info);
    if (ppc->ops->setfromoptions) {
      (*ppc->ops->setfromoptions)(ppc);
    }
    info = PCShellSetName(ppc, "bfgs"); CHKERRQ(info);
    info = PCShellSetContext(ppc, M); CHKERRQ(info);
    info = PCShellSetApply(ppc, bfgs_apply); CHKERRQ(info);
    break;

  default:
    // Use the pc method set by pc_type
    break;
  }

  // Initialize trust-region radius
  switch(tr->init_type) {
  case NTR_INIT_CONSTANT:
    // Use the initial radius specified
    break;

  case NTR_INIT_INTERPOLATION:
    // Use the initial radius specified
    max_radius = 0.0;

    for (j = 0; j < j_max; ++j) {
      fmin = f;
      sigma = 0.0;

      if (needH) {
        info = TaoComputeHessian(tao, X, H); CHKERRQ(info);
        needH = 0;
      }

      for (i = 0; i < i_max; ++i) {
	info = W->Waxpby(1.0, X, -radius / gnorm, G); CHKERRQ(info);

        info = TaoComputeFunction(tao, W, &ftrial); CHKERRQ(info);
        if (TaoInfOrNaN(ftrial)) {
	  tau = tr->gamma1_i;
        }
        else {
	  if (ftrial < fmin) {
            fmin = ftrial;
            sigma = -radius / gnorm;
          }

          info = H->Multiply(G, D); CHKERRQ(info);
          info = D->Dot(G, &prered); CHKERRQ(info);

          prered = radius * (gnorm - 0.5 * radius * prered / (gnorm * gnorm));
          actred = f - ftrial;
	  if ((fabs(actred) <= tr->epsilon) && 
              (fabs(prered) <= tr->epsilon)) {
	    kappa = 1.0;
	  }
	  else {
	    kappa = actred / prered;
	  }

	  tau_1 = tr->theta_i * gnorm * radius / (tr->theta_i * gnorm * radius + (1.0 - tr->theta_i) * prered - actred);
          tau_2 = tr->theta_i * gnorm * radius / (tr->theta_i * gnorm * radius - (1.0 + tr->theta_i) * prered + actred);
	  tau_min = TaoMin(tau_1, tau_2);
	  tau_max = TaoMax(tau_1, tau_2);

	  if (fabs(kappa - 1.0) <= tr->mu1_i) {
	    // Great agreement
            max_radius = TaoMax(max_radius, radius);

	    if (tau_max < 1.0) {
              tau = tr->gamma3_i;
            }
            else if (tau_max > tr->gamma4_i) {
              tau = tr->gamma4_i;
            }
            else {
              tau = tau_max;
            }
          }
          else if (fabs(kappa - 1.0) <= tr->mu2_i) {
	    // Good agreement
            max_radius = TaoMax(max_radius, radius);

	    if (tau_max < tr->gamma2_i) {
              tau = tr->gamma2_i;
            }
            else if (tau_max > tr->gamma3_i) {
              tau = tr->gamma3_i;
            }
            else {
              tau = tau_max;
            }
          }
          else {
	    // Not good agreement
	    if (tau_min > 1.0) {
	      tau = tr->gamma2_i;
            }
            else if (tau_max < tr->gamma1_i) {
              tau = tr->gamma1_i;
            }
	    else if ((tau_min < tr->gamma1_i) && (tau_max >= 1.0)) {
	      tau = tr->gamma1_i;
            }
	    else if ((tau_1 >= tr->gamma1_i) && (tau_1 < 1.0) && 
                     ((tau_2 < tr->gamma1_i) || (tau_2 >= 1.0))) {
              tau = tau_1;
            }
	    else if ((tau_2 >= tr->gamma1_i) && (tau_2 < 1.0) && 
                     ((tau_1 < tr->gamma1_i) || (tau_2 >= 1.0))) {
              tau = tau_2;
            }
            else {
              tau = tau_max;
            }
          }
        }
        radius = tau * radius;
      }
  
      if (fmin < f) {
        f = fmin;
        info = X->Axpy(sigma, G); CHKERRQ(info);
        info = TaoComputeGradient(tao, X, G); CHKERRQ(info);

        info = G->Norm2(&gnorm); CHKERRQ(info);
        if (TaoInfOrNaN(f) || TaoInfOrNaN(gnorm)) {
          SETERRQ(1, "User provided compute function generated Inf or NaN");
        }
        needH = 1;

        info = TaoMonitor(tao, iter, f, gnorm, 0.0, 1.0, &reason); CHKERRQ(info);
        if (reason != TAO_CONTINUE_ITERATING) {
          TaoFunctionReturn(0);
        }
      }
    }
    radius = TaoMax(radius, max_radius);

    // Modify the radius if it is too large or small
    radius = TaoMax(radius, tr->min_radius);
    radius = TaoMin(radius, tr->max_radius);
    break;

  default:
    // Norm of the first direction will initialize radius
    radius = 0.0;
    break;
  }

  // Set initial scaling for the BFGS preconditioner 
  // This step is done after computing the initial trust-region radius
  // since the function value may have decreased
  if (NTR_PC_BFGS == tr->pc_type) {
    if (f != 0.0) {
      delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
    }
    else {
      delta = 2.0 / (gnorm*gnorm);
    }
    info = M->SetDelta(delta); CHKERRQ(info);
  }

  // Have not converged; continue with Newton method
  while (reason == TAO_CONTINUE_ITERATING) {
    ++iter;

    // Compute the Hessian
    if (needH) {
      info = TaoComputeHessian(tao, X, H); CHKERRQ(info);
      needH = 0;
    }

    if (NTR_PC_BFGS == tr->pc_type) {
      if (BFGS_SCALE_AHESS == tr->bfgs_scale_type) {
        // Obtain diagonal for the bfgs preconditioner
        info = H->GetDiagonal(Diag); CHKERRQ(info);
        info = Diag->AbsoluteValue(); CHKERRQ(info);
        info = Diag->Reciprocal(); CHKERRQ(info);
        info = M->SetScale(Diag); CHKERRQ(info);
      }

      // Update the limited memory preconditioner
      info = M->Update(X, G); CHKERRQ(info);
      ++bfgsUpdates;
    }

    while (reason == TAO_CONTINUE_ITERATING) {
      // Solve the trust region subproblem
      info = TaoPreLinearSolve(tao, H); CHKERRQ(info);
      info = TaoLinearSolveTrustRegion(tao, H, G, D, radius, &success); CHKERRQ(info);
      info = pls->GetNormDirection(&norm_d); CHKERRQ(info);
      if (0.0 == radius) {
        // Radius was uninitialized; use the norm of the direction
        if (norm_d > 0.0) {
          radius = norm_d;

          // Modify the radius if it is too large or small
          radius = TaoMax(radius, tr->min_radius);
          radius = TaoMin(radius, tr->max_radius);
        }
        else {
          // The direction was bad; set radius to default value and re-solve 
          // the trust-region subproblem to get a direction
          info = TaoGetInitialTrustRegionRadius(tao, &radius); CHKERRQ(info);

          // Modify the radius if it is too large or small
          radius = TaoMax(radius, tr->min_radius);
          radius = TaoMin(radius, tr->max_radius);

          info = TaoLinearSolveTrustRegion(tao, H, G, D, radius, &success); CHKERRQ(info);
          info = pls->GetNormDirection(&norm_d); CHKERRQ(info);
	  if (norm_d == 0.0) {
            SETERRQ(1, "Initial direction zero");
          }
        }
      }
      info = D->Negate(); CHKERRQ(info);

      info = KSPGetConvergedReason(pksp, &ksp_reason); CHKERRQ(info);
      if ((KSP_DIVERGED_INDEFINITE_PC == ksp_reason) &&
          (NTR_PC_BFGS == tr->pc_type) && (bfgsUpdates > 1)) {
        // Preconditioner is numerically indefinite; reset the
        // approximate if using BFGS preconditioning.
  
        if (f != 0.0) {
          delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
        }
        else {
          delta = 2.0 / (gnorm*gnorm);
        }
        info = M->SetDelta(delta); CHKERRQ(info);
        info = M->Reset(); CHKERRQ(info);
        info = M->Update(X, G); CHKERRQ(info);
        bfgsUpdates = 1;
      }

      if (NTR_UPDATE_REDUCTION == tr->update_type) {
	// Get predicted reduction
	info = pls->GetObjFcn(&prered); CHKERRQ(info);
	
	if (prered >= 0.0) {
	  // The predicted reduction has the wrong sign.  This cannot
	  // happen in infinite precision arithmetic.  Step should
	  // be rejected!
	  radius = tr->alpha1 * TaoMin(radius, norm_d);
	}
	else {
	  // Compute trial step and function value
	  info = W->Waxpby(1.0, X, 1.0, D); CHKERRQ(info);
	  info = TaoComputeFunction(tao, W, &ftrial); CHKERRQ(info);
	  if (TaoInfOrNaN(ftrial)) {
	    radius = tr->alpha1 * TaoMin(radius, norm_d);
	  }
	  else {
	    // Compute and actual reduction
	    actred = f - ftrial;
	    prered = -prered;
	    if ((fabs(actred) <= tr->epsilon) && 
                (fabs(prered) <= tr->epsilon)) {
	      kappa = 1.0;
	    }
	    else {
	      kappa = actred / prered;
	    }
	    
	    // Accept or reject the step and update radius
	    if (kappa < tr->eta1) {
	      // Reject the step
	      radius = tr->alpha1 * TaoMin(radius, norm_d);
	    } 
	    else {
	      // Accept the step
	      if (kappa < tr->eta2) { 
		// Marginal bad step
		radius = tr->alpha2 * TaoMin(radius, norm_d);
	      }
	      else if (kappa < tr->eta3) {
		// Reasonable step
		radius = tr->alpha3 * radius;
	      }
	      else if (kappa < tr->eta4) { 
		// Good step
		radius = TaoMax(tr->alpha4 * norm_d, radius);
	      }
	      else {
		// Very good step
		radius = TaoMax(tr->alpha5 * norm_d, radius);
	      }
	      break;
	    }
	  }
	} 
      }
      else {
	// Get predicted reduction
	info = pls->GetObjFcn(&prered); CHKERRQ(info);

	if (prered >= 0.0) {
	  // The predicted reduction has the wrong sign.  This cannot
	  // happen in infinite precision arithmetic.  Step should
	  // be rejected!
	  radius = tr->gamma1 * TaoMin(radius, norm_d);
	}
	else {
	  info = W->Waxpby(1.0, X, 1.0, D); CHKERRQ(info);
	  info = TaoComputeFunction(tao, W, &ftrial); CHKERRQ(info);
	  if (TaoInfOrNaN(ftrial)) {
	    radius = tr->gamma1 * TaoMin(radius, norm_d);
	  }
	  else {
	    info = D->Dot(G, &beta); CHKERRQ(info);

	    actred = f - ftrial;
	    prered = -prered;
	    if ((fabs(actred) <= tr->epsilon) && 
                (fabs(prered) <= tr->epsilon)) {
	      kappa = 1.0;
	    }
	    else {
	      kappa = actred / prered;
	    }

	    tau_1 = tr->theta * beta / (tr->theta * beta - (1.0 - tr->theta) * prered + actred);
	    tau_2 = tr->theta * beta / (tr->theta * beta + (1.0 + tr->theta) * prered - actred);
	    tau_min = TaoMin(tau_1, tau_2);
	    tau_max = TaoMax(tau_1, tau_2);

	    if (kappa >= 1.0 - tr->mu1) {
	      // Great agreement; accept step and update radius
	      if (tau_max < 1.0) {
		radius = TaoMax(radius, tr->gamma3 * norm_d);
	      }
	      else if (tau_max > tr->gamma4) {
		radius = TaoMax(radius, tr->gamma4 * norm_d);
	      }
	      else {
		radius = TaoMax(radius, tau_max * norm_d);
	      }
	      break;
	    }
	    else if (kappa >= 1.0 - tr->mu2) {
	      // Good agreement

	      if (tau_max < tr->gamma2) {
		radius = tr->gamma2 * TaoMin(radius, norm_d);
	      }
	      else if (tau_max > tr->gamma3) {
		radius = TaoMax(radius, tr->gamma3 * norm_d);
	      }
	      else if (tau_max < 1.0) {
		radius = tau_max * TaoMin(radius, norm_d);
	      }
	      else {
		radius = TaoMax(radius, tau_max * norm_d);
	      }
	      break;
	    }
	    else {
	      // Not good agreement
	      if (tau_min > 1.0) {
		radius = tr->gamma2 * TaoMin(radius, norm_d);
	      }
	      else if (tau_max < tr->gamma1) {
		radius = tr->gamma1 * TaoMin(radius, norm_d);
	      }
	      else if ((tau_min < tr->gamma1) && (tau_max >= 1.0)) {
		radius = tr->gamma1 * TaoMin(radius, norm_d);
	      }
	      else if ((tau_1 >= tr->gamma1) && (tau_1 < 1.0) && 
		       ((tau_2 < tr->gamma1) || (tau_2 >= 1.0))) {
		radius = tau_1 * TaoMin(radius, norm_d);
	      }
	      else if ((tau_2 >= tr->gamma1) && (tau_2 < 1.0) && 
		       ((tau_1 < tr->gamma1) || (tau_2 >= 1.0))) {
		radius = tau_2 * TaoMin(radius, norm_d);
	      }
	      else {
		radius = tau_max * TaoMin(radius, norm_d);
	      }
	    }
	  }
	}
      }

      // The step computed was not good and the radius was decreased.
      // Monitor the radius to terminate.
      info = TaoMonitor(tao, iter, f, gnorm, 0.0, radius, &reason); CHKERRQ(info);
    }

    // The radius may have been increased; modify if it is too large
    radius = TaoMin(radius, tr->max_radius);

    if (reason == TAO_CONTINUE_ITERATING) {
      info = X->CopyFrom(W); CHKERRQ(info);
      f = ftrial;
      info = TaoComputeGradient(tao, X, G); CHKERRQ(info);
      info = G->Norm2(&gnorm); CHKERRQ(info);
      if (TaoInfOrNaN(f) || TaoInfOrNaN(gnorm)) {
	SETERRQ(1, "User provided compute function generated Inf or NaN");
      }
      needH = 1;

      info = TaoMonitor(tao, iter, f, gnorm, 0.0, radius, &reason); CHKERRQ(info);
    }
  }
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_NTR"
static int TaoSetUp_NTR(TAO_SOLVER tao, void *solver)
{
  TAO_NTR *tr = (TAO_NTR *)solver;
  TaoVec *X;
  TaoMat *H;
  int info;

  TaoFunctionBegin;

  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = X->Clone(&tr->G); CHKERRQ(info);
  info = X->Clone(&tr->D); CHKERRQ(info);
  info = X->Clone(&tr->W); CHKERRQ(info);

  tr->Diag = 0;
  tr->M = 0;

  info = TaoSetLagrangianGradientVector(tao, tr->G); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, tr->D); CHKERRQ(info);

  // Set linear solver to default for trust region
  info = TaoGetHessian(tao, &H); CHKERRQ(info);
  info = TaoCreateLinearSolver(tao, H, 200, 0); CHKERRQ(info); 

  // Check sizes for compatability
  info = TaoCheckFGH(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_NTR"
static int TaoDestroy_NTR(TAO_SOLVER tao, void *solver)
{
  TAO_NTR *tr = (TAO_NTR *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoVecDestroy(tr->G); CHKERRQ(info);
  info = TaoVecDestroy(tr->D); CHKERRQ(info);
  info = TaoVecDestroy(tr->W); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao, 0); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, 0); CHKERRQ(info);

  if (tr->Diag) {
    info = TaoVecDestroy(tr->Diag); CHKERRQ(info);
    tr->Diag = 0;
  }

  if (tr->M) {
    info = TaoMatDestroy(tr->M); CHKERRQ(info);
    tr->M = 0;
  }
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_NTR"
static int TaoSetOptions_NTR(TAO_SOLVER tao, void*solver)
{
  TAO_NTR *tr = (TAO_NTR *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoOptionsHead("Newton trust region method for unconstrained optimization"); CHKERRQ(info);
  info = TaoOptionList("-tao_ntr_ksp_type", "ksp type", "", NTR_KSP, NTR_KSP_TYPES, NTR_KSP[tr->ksp_type], &tr->ksp_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_ntr_pc_type", "pc type", "", NTR_PC, NTR_PC_TYPES, NTR_PC[tr->pc_type], &tr->pc_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_ntr_bfgs_scale_type", "bfgs scale type", "", BFGS_SCALE, BFGS_SCALE_TYPES, BFGS_SCALE[tr->bfgs_scale_type], &tr->bfgs_scale_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_ntr_init_type", "radius initialization type", "", NTR_INIT, NTR_INIT_TYPES, NTR_INIT[tr->init_type], &tr->init_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_ntr_update_type", "radius update type", "", NTR_UPDATE, NTR_UPDATE_TYPES, NTR_UPDATE[tr->update_type], &tr->update_type, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_eta1", "step is unsuccessful if actual reduction < eta1 * predicted reduction", "", tr->eta1, &tr->eta1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_eta2", "", "", tr->eta2, &tr->eta2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_eta3", "", "", tr->eta3, &tr->eta3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_eta4", "", "", tr->eta4, &tr->eta4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_alpha1", "", "", tr->alpha1, &tr->alpha1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_alpha2", "", "", tr->alpha2, &tr->alpha2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_alpha3", "", "", tr->alpha3, &tr->alpha3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_alpha4", "", "", tr->alpha4, &tr->alpha4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_alpha5", "", "", tr->alpha5, &tr->alpha5, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_mu1", "", "", tr->mu1, &tr->mu1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_mu2", "", "", tr->mu2, &tr->mu2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_gamma1", "", "", tr->gamma1, &tr->gamma1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_gamma2", "", "", tr->gamma2, &tr->gamma2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_gamma3", "", "", tr->gamma3, &tr->gamma3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_gamma4", "", "", tr->gamma4, &tr->gamma4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_theta", "", "", tr->theta, &tr->theta, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_mu1_i", "", "", tr->mu1_i, &tr->mu1_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_mu2_i", "", "", tr->mu2_i, &tr->mu2_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_gamma1_i", "", "", tr->gamma1_i, &tr->gamma1_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_gamma2_i", "", "", tr->gamma2_i, &tr->gamma2_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_gamma3_i", "", "", tr->gamma3_i, &tr->gamma3_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_gamma4_i", "", "", tr->gamma4_i, &tr->gamma4_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_theta_i", "", "", tr->theta_i, &tr->theta_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_min_radius", "lower bound on initial trust-region radius", "", tr->min_radius, &tr->min_radius, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_max_radius", "upper bound on trust-region radius", "", tr->max_radius, &tr->max_radius, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntr_epsilon", "tolerance used when computing actual and predicted reduction", "", tr->epsilon, &tr->epsilon, 0); CHKERRQ(info);
  info = TaoOptionsTail(); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_NTR"
static int TaoView_NTR(TAO_SOLVER tao,void*solver)
{
  TAO_NTR *tr = (TAO_NTR *)solver;
  int info;

  TaoFunctionBegin;
  if (NTR_PC_BFGS == tr->pc_type && tr->M) {
    info = TaoPrintInt(tao, "  Rejected matrix updates: %d\n", tr->M->GetRejects()); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_NTR"
int TaoCreate_NTR(TAO_SOLVER tao)
{
  TAO_NTR *tr;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_NTR, &tr); CHKERRQ(info);
  info = PetscLogObjectMemory(tao, sizeof(TAO_NTR)); CHKERRQ(info);

  info = TaoSetTaoSolveRoutine(tao, TaoSolve_NTR, (void *)tr); CHKERRQ(info);
  info = TaoSetTaoSetUpDownRoutines(tao, TaoSetUp_NTR, TaoDestroy_NTR); CHKERRQ(info);
  info = TaoSetTaoOptionsRoutine(tao, TaoSetOptions_NTR); CHKERRQ(info);
  info = TaoSetTaoViewRoutine(tao, TaoView_NTR); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao, 50); CHKERRQ(info);
  info = TaoSetTolerances(tao, 1e-10, 1e-10, 0, 0); CHKERRQ(info);

  info = TaoSetTrustRegionRadius(tao, 100.0); CHKERRQ(info);
  info = TaoSetTrustRegionTolerance(tao, 1.0e-12); CHKERRQ(info);

  // Standard trust region update parameters
  tr->eta1 = 1.0e-4;
  tr->eta2 = 0.25;
  tr->eta3 = 0.50;
  tr->eta4 = 0.90;

  tr->alpha1 = 0.25;
  tr->alpha2 = 0.50;
  tr->alpha3 = 1.00;
  tr->alpha4 = 2.00;
  tr->alpha5 = 4.00;

  // Interpolation parameters
  tr->mu1_i = 0.35;
  tr->mu2_i = 0.50;

  tr->gamma1_i = 0.0625;
  tr->gamma2_i = 0.50;
  tr->gamma3_i = 2.00;
  tr->gamma4_i = 5.00;

  tr->theta_i = 0.25;

  // Interpolation trust region update parameters
  tr->mu1 = 0.10;
  tr->mu2 = 0.50;

  tr->gamma1 = 0.25;
  tr->gamma2 = 0.50;
  tr->gamma3 = 2.00;
  tr->gamma4 = 4.00;

  tr->theta = 0.05;

  tr->min_radius = 1.0e-10;
  tr->max_radius = 1.0e10;
  tr->epsilon = 1.0e-6;

  tr->ksp_type        = NTR_KSP_STCG;
  tr->pc_type         = NTR_PC_BFGS;
  tr->bfgs_scale_type = BFGS_SCALE_AHESS;
  tr->init_type	      = NTR_INIT_INTERPOLATION;
  tr->update_type     = NTR_UPDATE_REDUCTION;
  TaoFunctionReturn(0);
}
EXTERN_C_END

#endif

