/*$Id: ntl.c,v 1.4 2009-08-25 16:25:11 sarich Exp $*/

#include "ntl.h"             /*I "tao_solver.h" I*/

#ifdef TAO_USE_PETSC
#include "petscksp.h"
#include "petscpc.h"
#include "src/petsctao/linearsolver/taolinearsolver_petsc.h"
#include "src/petsctao/vector/taovec_petsc.h"

#include "private/kspimpl.h"
#include "private/pcimpl.h"

#define NTL_KSP_NASH	0
#define NTL_KSP_STCG	1
#define NTL_KSP_GLTR	2
#define NTL_KSP_TYPES	3

#define NTL_PC_NONE	0
#define NTL_PC_AHESS	1
#define NTL_PC_BFGS	2
#define NTL_PC_PETSC	3
#define NTL_PC_TYPES	4

#define BFGS_SCALE_AHESS	0
#define BFGS_SCALE_BFGS		1
#define BFGS_SCALE_TYPES	2

#define NTL_INIT_CONSTANT         0
#define NTL_INIT_DIRECTION        1
#define NTL_INIT_INTERPOLATION    2
#define NTL_INIT_TYPES            3

#define NTL_UPDATE_REDUCTION      0
#define NTL_UPDATE_INTERPOLATION  1
#define NTL_UPDATE_TYPES          2

static const char *NTL_KSP[64] = {
  "nash", "stcg", "gltr"
};

static const char *NTL_PC[64] = {
  "none", "ahess", "bfgs", "petsc"
};

static const char *BFGS_SCALE[64] = {
  "ahess", "bfgs"
};

static const char *NTL_INIT[64] = {
  "constant", "direction", "interpolation"
};

static const char *NTL_UPDATE[64] = {
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

// Implements Newton's Method with a trust-region, line-search approach for 
// solving unconstrained minimization problems.  A More'-Thuente line search 
// is used to guarantee that the bfgs preconditioner remains positive
// definite.

#define NTL_NEWTON 		0
#define NTL_BFGS 		1
#define NTL_SCALED_GRADIENT 	2
#define NTL_GRADIENT 		3

#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_NTL"
static int TaoSolve_NTL(TAO_SOLVER tao, void *solver)
{
  TAO_NTL *tl = (TAO_NTL *)solver;
  TaoVec *X, *G = tl->G, *D = tl->D, *W = tl->W, *PG=tl->PG;
  TaoVec *Xold = tl->Xold, *Gold = tl->Gold, *Diag = tl->Diag;
  TaoMat *H;
  TaoLMVMMat *M = tl->M;

  TaoLinearSolver *tls;
  TaoLinearSolverPetsc *pls;

  KSP pksp;
  PC ppc;

  KSPConvergedReason ksp_reason;
  TaoTerminateReason reason;
  TaoTruth success;
  
  double fmin, ftrial, f_full, prered, actred, kappa, sigma;
  double tau, tau_1, tau_2, tau_max, tau_min, max_radius;
  double f, fold, gdx, gnorm;
  double step = 1.0;

  double delta;
  double radius, norm_d = 0.0;

  int info;
  TaoInt stepType;
  TaoInt iter = 0, status = 0;
  TaoInt bfgsUpdates = 0;
  TaoInt needH;

  TaoInt i_max = 5;
  TaoInt j_max = 1;
  TaoInt i, j;

  TaoInt tr_reject;

  TaoFunctionBegin;

  // Initialize trust-region radius
  info = TaoGetInitialTrustRegionRadius(tao, &radius); CHKERRQ(info);
  if (radius < 0.0) {
    SETERRQ(1, "Initial radius negative");
  }

  // Modify the radius if it is too large or small
  radius = TaoMax(radius, tl->min_radius);
  radius = TaoMin(radius, tl->max_radius);

  // Get vectors we will need
  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = TaoGetHessian(tao, &H); CHKERRQ(info);

  /*   Project the current point onto the feasible set */
  TaoVec *XU, *XL;
  info = TaoGetVariableBounds(tao,&XL,&XU);CHKERRQ(info);
  info = TaoEvaluateVariableBounds(tao,XL,XU); CHKERRQ(info);
  info = X->Median(XL,X,XU); CHKERRQ(info);
  
  if (NTL_PC_BFGS == tl->pc_type && !M) {
    tl->M = new TaoLMVMMat(X);
    M = tl->M;
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
  if ((NTL_PC_BFGS == tl->pc_type) && 
      (BFGS_SCALE_BFGS != tl->bfgs_scale_type)) {
    if (!Diag) {
      info = X->Clone(&tl->Diag); CHKERRQ(info);
      Diag = tl->Diag;
    }
  }

  // Modify the linear solver to a conjugate gradient method
  info = TaoGetLinearSolver(tao, &tls); CHKERRQ(info);
  pls  = dynamic_cast <TaoLinearSolverPetsc *> (tls);

  pksp = pls->GetKSP();
  switch(tl->ksp_type) {
  case NTL_KSP_NASH:
    info = KSPSetType(pksp, KSPNASH); CHKERRQ(info);
    if (pksp->ops->setfromoptions) {
      (*pksp->ops->setfromoptions)(pksp);
    }
    break;

  case NTL_KSP_STCG:
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
  switch(tl->pc_type) {
  case NTL_PC_NONE:
    info = PCSetType(ppc, PCNONE); CHKERRQ(info);
    if (ppc->ops->setfromoptions) {
      (*ppc->ops->setfromoptions)(ppc);
    }
    break;

  case NTL_PC_AHESS:
    info = PCSetType(ppc, PCJACOBI); CHKERRQ(info);
    if (ppc->ops->setfromoptions) {
      (*ppc->ops->setfromoptions)(ppc);
    }
    info = PCJacobiSetUseAbs(ppc); CHKERRQ(info);
    break;

  case NTL_PC_BFGS:
    info = KSPSetNormType(pksp, KSP_NORM_PRECONDITIONED); CHKERRQ(info);
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

  // Initialize trust-region radius.  The initialization is only performed 
  // when we are using Steihaug-Toint or the Generalized Lanczos method.
  info=PetscInfo1(tao,"initializing trust region radius %d \n",tl->init_type);
  switch(tl->init_type) {
  case NTL_INIT_CONSTANT:
    // Use the initial radius specified
    break;

  case NTL_INIT_INTERPOLATION:
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
          tau = tl->gamma1_i;
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
          if ((fabs(actred) <= tl->epsilon) && 
              (fabs(prered) <= tl->epsilon)) {
            kappa = 1.0;
          }
          else {
            kappa = actred / prered;
          }

          tau_1 = tl->theta_i * gnorm * radius / (tl->theta_i * gnorm * radius + (1.0 - tl->theta_i) * prered - actred);
          tau_2 = tl->theta_i * gnorm * radius / (tl->theta_i * gnorm * radius - (1.0 + tl->theta_i) * prered + actred);
          tau_min = TaoMin(tau_1, tau_2);
          tau_max = TaoMax(tau_1, tau_2);

          if (fabs(kappa - 1.0) <= tl->mu1_i) {
            // Great agreement
            max_radius = TaoMax(max_radius, radius);

            if (tau_max < 1.0) {
              tau = tl->gamma3_i;
            }
            else if (tau_max > tl->gamma4_i) {
              tau = tl->gamma4_i;
            }
            else if (tau_1 >= 1.0 && tau_1 <= tl->gamma4_i && tau_2 < 1.0) {
              tau = tau_1;
            }
            else if (tau_2 >= 1.0 && tau_2 <= tl->gamma4_i && tau_1 < 1.0) {
              tau = tau_2;
            }
            else {
              tau = tau_max;
            }
          }
          else if (fabs(kappa - 1.0) <= tl->mu2_i) {
            // Good agreement
            max_radius = TaoMax(max_radius, radius);

            if (tau_max < tl->gamma2_i) {
	      tau = tl->gamma2_i;
	    }
	    else if (tau_max > tl->gamma3_i) {
	      tau = tl->gamma3_i;
	    }
	    else {
	      tau = tau_max;
	    }
	  }
	  else {
	    // Not good agreement
	    if (tau_min > 1.0) {
	      tau = tl->gamma2_i;
	    }
	    else if (tau_max < tl->gamma1_i) {
	      tau = tl->gamma1_i;
	    }
	    else if ((tau_min < tl->gamma1_i) && (tau_max >= 1.0)) {
	      tau = tl->gamma1_i;
	    }
	    else if ((tau_1 >= tl->gamma1_i) && (tau_1 < 1.0) &&
		     ((tau_2 < tl->gamma1_i) || (tau_2 >= 1.0))) {
	      tau = tau_1;
	    }
	    else if ((tau_2 >= tl->gamma1_i) && (tau_2 < 1.0) &&
		     ((tau_1 < tl->gamma1_i) || (tau_2 >= 1.0))) {
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
    radius = TaoMax(radius, tl->min_radius);
    radius = TaoMin(radius, tl->max_radius);
    break;

  default:
    // Norm of the first direction will initialize radius
    radius = 0.0;
    break;
  }
  info=PetscInfo1(tao,"initial trust region radius %22.15e \n",radius);

  // Set initial scaling for the BFGS preconditioner
  // This step is done after computing the initial trust-region radius
  // since the function value may have decreased
  if (NTL_PC_BFGS == tl->pc_type) {
    if (f != 0.0) {
      delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
    }
    else {
      delta = 2.0 / (gnorm*gnorm);
    }
    info = M->SetDelta(delta); CHKERRQ(info);
  }

  // Set counter for gradient/reset steps
  tl->trust = 0;
  tl->newt = 0;
  tl->bfgs = 0;
  tl->sgrad = 0;
  tl->grad = 0;

  // Have not converged; continue with Newton method
  while (reason == TAO_CONTINUE_ITERATING) {
    ++iter;

    // Compute the Hessian
    if (needH) {
      info = TaoComputeHessian(tao, X, H); CHKERRQ(info);
      needH = 0;
    }

    if (NTL_PC_BFGS == tl->pc_type) {
      if (BFGS_SCALE_AHESS == tl->bfgs_scale_type) {
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

    // use constrained norm for tolerance
    PetscScalar boundNorm;
    info = PG->BoundGradientProjection(G,XL,X,XU);CHKERRQ(info);
    info = PG->Norm2(&boundNorm); CHKERRQ(info);
    
    PetscScalar ewAtol  = PetscMin(0.5,boundNorm)*boundNorm;
    info = KSPSetTolerances(pksp,PETSC_DEFAULT,ewAtol,
                            PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(info);
    info=PetscInfo2(tao,"TaoSolve_NTL: gnorm=%22.12e, boundNorm=%22.12e\n",gnorm,boundNorm);
    pksp->printreason = PETSC_TRUE;
    info = KSPView(pksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(info);

    // Solve the Newton system of equations
    info = TaoPreLinearSolve(tao, H); CHKERRQ(info);
    info = TaoLinearSolveTrustRegion(tao, H, G, D, radius, &success); CHKERRQ(info);
    info = pls->GetNormDirection(&norm_d); CHKERRQ(info);
    if (0.0 == radius) {
      // Radius was uninitialized; use the norm of the direction
      if (norm_d > 0.0) {
	radius = norm_d;

	// Modify the radius if it is too large or small
	radius = TaoMax(radius, tl->min_radius);
	radius = TaoMin(radius, tl->max_radius);
      }
      else {
	// The direction was bad; set radius to default value and re-solve 
	// the trust-region subproblem to get a direction
	info = TaoGetInitialTrustRegionRadius(tao, &radius); CHKERRQ(info);

	// Modify the radius if it is too large or small
	radius = TaoMax(radius, tl->min_radius);
	radius = TaoMin(radius, tl->max_radius);

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
        (NTL_PC_BFGS == tl->pc_type) && (bfgsUpdates > 1)) {
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

    // Check trust-region reduction conditions
    tr_reject = 0;
    info=PetscInfo1(tao,"ntl update type , %d\n",tl->update_type);
    if (NTL_UPDATE_REDUCTION == tl->update_type) {
      // Get predicted reduction
      info = pls->GetObjFcn(&prered); CHKERRQ(info);
      if (prered >= 0.0) {
	// The predicted reduction has the wrong sign.  This cannot
	// happen in infinite precision arithmetic.  Step should
	// be rejected!
	radius = tl->alpha1 * TaoMin(radius, norm_d);
	tr_reject = 1;
      }
      else {
	// Compute trial step and function value
	info = W->Waxpby(1.0, X, 1.0, D); CHKERRQ(info);
	info = TaoComputeFunction(tao, W, &ftrial); CHKERRQ(info);
	if (TaoInfOrNaN(ftrial)) {
	  radius = tl->alpha1 * TaoMin(radius, norm_d);
	  tr_reject = 1;
	}
	else {
	  // Compute and actual reduction
	  actred = f - ftrial;
	  prered = -prered;
	  if ((fabs(actred) <= tl->epsilon) &&
	      (fabs(prered) <= tl->epsilon)) {
	    kappa = 1.0;
	  }
	  else {
	    kappa = actred / prered;
	  }
          info=PetscInfo1(tao,"Checking TR model kappa %22.15e \n",kappa );

	  // Accept of reject the step and update radius
	  if (kappa < tl->eta1) {
	    // Reject the step
	    radius = tl->alpha1 * TaoMin(radius, norm_d);
	    tr_reject = 1;
            info=PetscInfo1(tao,"Reject the TR step radius %22.15e \n",radius);
	  }
	  else {
	    // Accept the step
	    if (kappa < tl->eta2) {
	      // Marginal bad step
	      radius = tl->alpha2 * TaoMin(radius, norm_d);
              info=PetscInfo1(tao,"Marginal bad TR step %22.15e \n",radius);
	    }
	    else if (kappa < tl->eta3) {
	      // Reasonable step
	      radius = tl->alpha3 * radius;
              info=PetscInfo1(tao,"Reasonable TR step %22.15e \n",radius);
	    }
	    else if (kappa < tl->eta4) {
	      // Good step
	      radius = TaoMax(tl->alpha4 * norm_d, radius);
              info=PetscInfo1(tao,"Good TR step %22.15e \n",radius);
	    }
	    else {
	      // Very good step
	      radius = TaoMax(tl->alpha5 * norm_d, radius);
              info=PetscInfo1(tao,"Very good TR step %22.15e \n",radius);
	    }
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
	radius = tl->gamma1 * TaoMin(radius, norm_d);
	tr_reject = 1;
      }
      else {
	info = W->Waxpby(1.0, X, 1.0, D); CHKERRQ(info);
	info = TaoComputeFunction(tao, W, &ftrial); CHKERRQ(info);
	if (TaoInfOrNaN(ftrial)) {
	  radius = tl->gamma1 * TaoMin(radius, norm_d);
	  tr_reject = 1;
	}
	else {
	  info = D->Dot(G, &gdx); CHKERRQ(info);

	  actred = f - ftrial;
	  prered = -prered;
	  if ((fabs(actred) <= tl->epsilon) &&
	      (fabs(prered) <= tl->epsilon)) {
	    kappa = 1.0;
	  }
	  else {
	    kappa = actred / prered;
	  }

	  tau_1 = tl->theta * gdx / (tl->theta * gdx - (1.0 - tl->theta) * prered + actred);
	  tau_2 = tl->theta * gdx / (tl->theta * gdx + (1.0 + tl->theta) * prered - actred);
	  tau_min = TaoMin(tau_1, tau_2);
	  tau_max = TaoMax(tau_1, tau_2);

	  if (kappa >= 1.0 - tl->mu1) {
	    // Great agreement; accept step and update radius
	    if (tau_max < 1.0) {
	      radius = TaoMax(radius, tl->gamma3 * norm_d);
	    }
	    else if (tau_max > tl->gamma4) {
	      radius = TaoMax(radius, tl->gamma4 * norm_d);
	    }
	    else {
	      radius = TaoMax(radius, tau_max * norm_d);
	    }
	  }
	  else if (kappa >= 1.0 - tl->mu2) {
	    // Good agreement

	    if (tau_max < tl->gamma2) {
	      radius = tl->gamma2 * TaoMin(radius, norm_d);
	    }
	    else if (tau_max > tl->gamma3) {
	      radius = TaoMax(radius, tl->gamma3 * norm_d);
	    }              else if (tau_max < 1.0) {
	      radius = tau_max * TaoMin(radius, norm_d);
	    }
	    else {
	      radius = TaoMax(radius, tau_max * norm_d);
	    }
	  }
	  else {
	    // Not good agreement
	    if (tau_min > 1.0) {
	      radius = tl->gamma2 * TaoMin(radius, norm_d);
	    }
	    else if (tau_max < tl->gamma1) {
	      radius = tl->gamma1 * TaoMin(radius, norm_d);
	    }
	    else if ((tau_min < tl->gamma1) && (tau_max >= 1.0)) {
	      radius = tl->gamma1 * TaoMin(radius, norm_d);
	    }
	    else if ((tau_1 >= tl->gamma1) && (tau_1 < 1.0) &&
		     ((tau_2 < tl->gamma1) || (tau_2 >= 1.0))) {
	      radius = tau_1 * TaoMin(radius, norm_d);
	    }
	    else if ((tau_2 >= tl->gamma1) && (tau_2 < 1.0) &&
		     ((tau_1 < tl->gamma1) || (tau_2 >= 1.0))) {
	      radius = tau_2 * TaoMin(radius, norm_d);
	    }
	    else {
	      radius = tau_max * TaoMin(radius, norm_d);
	    }
	    tr_reject = 1;
	  }
	}
      }
    }

    if (tr_reject) {
      // The trust-region constraints rejected the step.  Apply a linesearch.
      // Check for descent direction.
      info = D->Dot(G, &gdx); CHKERRQ(info);
      info=PetscInfo1(tao,"trust-region constraints rejected the step, gdx %22.15e \n",gdx);
      if ((gdx >= 0.0) || TaoInfOrNaN(gdx)) {
	// Newton step is not descent or direction produced Inf or NaN
	
	if (NTL_PC_BFGS != tl->pc_type) {
	  // We don't have the bfgs matrix around and updated
	  // Must use gradient direction in this case
	  info = D->CopyFrom(G); CHKERRQ(info);
	  info = D->Negate(); CHKERRQ(info);
	  
	  ++tl->grad;
	  stepType = NTL_GRADIENT;
	}
	else {
	  // Attempt to use the BFGS direction
	  info = M->Solve(G, D, &success); CHKERRQ(info);
	  info = D->Negate(); CHKERRQ(info);
	  
	  // Check for success (descent direction)
	  info = D->Dot(G, &gdx); CHKERRQ(info);
          info=PetscInfo1(tao,"newton direction fail bfgs descent gdx %22.15e \n",gdx);
	  if ((gdx >= 0) || TaoInfOrNaN(gdx)) {
	    // BFGS direction is not descent or direction produced not a number
	    // We can assert bfgsUpdates > 1 in this case because
	    // the first solve produces the scaled gradient direction,
	    // which is guaranteed to be descent
	    
	    // Use steepest descent direction (scaled)
	    
	    if (f != 0.0) {
	      delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
	    }
	    else {
	      delta = 2.0 / (gnorm*gnorm);
	    }
	    info = M->SetDelta(delta); CHKERRQ(info);
	    info = M->Reset(); CHKERRQ(info);
	    info = M->Update(X, G); CHKERRQ(info);
	    info = M->Solve(G, D, &success); CHKERRQ(info);
	    info = D->Negate(); CHKERRQ(info);
	    
	    bfgsUpdates = 1;
	    ++tl->sgrad;
	    stepType = NTL_SCALED_GRADIENT;
	  }
	  else {
	    if (1 == bfgsUpdates) {
	      // The first BFGS direction is always the scaled gradient
	      ++tl->sgrad;
	      stepType = NTL_SCALED_GRADIENT;
	    }
	    else {
	      ++tl->bfgs;
	      stepType = NTL_BFGS;
	    }
	  }
	}
      }
      else {
	// Computed Newton step is descent
        info=PetscInfo1(tao,"Computed Newton step is descent, gdx %22.15e \n",gdx);
	++tl->newt;
	stepType = NTL_NEWTON;
      }
      
      // Perform the linesearch
      fold = f;
      info = Xold->CopyFrom(X); CHKERRQ(info);
      info = Gold->CopyFrom(G); CHKERRQ(info);

      step = 1.0;
      info = TaoLineSearchApply(tao, X, G, D, W, &f, &f_full, &step, &status); CHKERRQ(info);

      while (status && stepType != NTL_GRADIENT) {
        info=PetscInfo1(tao,"line search fail , %d\n",status);
	// Linesearch failed
	f = fold;
	info = X->CopyFrom(Xold); CHKERRQ(info);
	info = G->CopyFrom(Gold); CHKERRQ(info);
	
	switch(stepType) {
	case NTL_NEWTON:
	  // Failed to obtain acceptable iterate with Newton step

	  if (NTL_PC_BFGS != tl->pc_type) {
	    // We don't have the bfgs matrix around and being updated
	    // Must use gradient direction in this case
	    info = D->CopyFrom(G); CHKERRQ(info);
	    
	    ++tl->grad;
	    stepType = NTL_GRADIENT;
	  }
	  else {
	    // Attempt to use the BFGS direction
	    info = M->Solve(G, D, &success); CHKERRQ(info);
	    
	    // Check for success (descent direction)
	    info = D->Dot(G, &gdx); CHKERRQ(info);
	    if ((gdx <= 0) || TaoInfOrNaN(gdx)) {
	      // BFGS direction is not descent or direction produced 
	      // not a number.  We can assert bfgsUpdates > 1 in this case
	      // Use steepest descent direction (scaled)
    
	      if (f != 0.0) {
		delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
	      }
	      else {
		delta = 2.0 / (gnorm*gnorm);
	      }
	      info = M->SetDelta(delta); CHKERRQ(info);
	      info = M->Reset(); CHKERRQ(info);
	      info = M->Update(X, G); CHKERRQ(info);
	      info = M->Solve(G, D, &success); CHKERRQ(info);
	      
	      bfgsUpdates = 1;
	      ++tl->sgrad;
	      stepType = NTL_SCALED_GRADIENT;
	    }
	    else {
	      if (1 == bfgsUpdates) {
		// The first BFGS direction is always the scaled gradient
		++tl->sgrad;
		stepType = NTL_SCALED_GRADIENT;
	      }
	      else {
		++tl->bfgs;
		stepType = NTL_BFGS;
	      }
	    }
	  }
	  break;

	case NTL_BFGS:
	  // Can only enter if pc_type == NTL_PC_BFGS
	  // Failed to obtain acceptable iterate with BFGS step
	  // Attempt to use the scaled gradient direction
	  
	  if (f != 0.0) {
	    delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
	  }
	  else {
	    delta = 2.0 / (gnorm*gnorm);
	  }
	  info = M->SetDelta(delta); CHKERRQ(info);
	  info = M->Reset(); CHKERRQ(info);
	  info = M->Update(X, G); CHKERRQ(info);
	  info = M->Solve(G, D, &success); CHKERRQ(info);
	  
	  bfgsUpdates = 1;
	  ++tl->sgrad;
	  stepType = NTL_SCALED_GRADIENT;
	  break;
	  
	case NTL_SCALED_GRADIENT:
	  // Can only enter if pc_type == NTL_PC_BFGS
	  // The scaled gradient step did not produce a new iterate;
	  // attemp to use the gradient direction.
	  // Need to make sure we are not using a different diagonal scaling
	  info = M->SetScale(0); CHKERRQ(info);
	  info = M->SetDelta(1.0); CHKERRQ(info);
	  info = M->Reset(); CHKERRQ(info);
	  info = M->Update(X, G); CHKERRQ(info);
	  info = M->Solve(G, D, &success); CHKERRQ(info);
	  
	  bfgsUpdates = 1;
	  ++tl->grad;
	  stepType = NTL_GRADIENT;
	  break;
	}
	info = D->Negate(); CHKERRQ(info);

	// This may be incorrect; linesearch has values for stepmax and stepmin
	// that should be reset.
	step = 1.0;
	info = TaoLineSearchApply(tao, X, G, D, W, &f, &f_full, &step, &status); CHKERRQ(info);
      }

      if (status) {
	// Failed to find an improving point
	f = fold;
	info = X->CopyFrom(Xold); CHKERRQ(info);
	info = G->CopyFrom(Gold); CHKERRQ(info);
	radius = 0.0;
      }
      else if (stepType == NTL_NEWTON) {
	if (step < tl->nu1) {
	  // Very bad step taken; reduce radius
	  radius = tl->omega1 * TaoMin(norm_d, radius);
          info=PetscInfo1(tao,"Very bad step taken; reduce radius %22.15e \n",radius);
	}
	else if (step < tl->nu2) {
	  // Reasonably bad step taken; reduce radius
	  radius = tl->omega2 * TaoMin(norm_d, radius);
          info=PetscInfo1(tao,"Reasonably bad step taken; reduce radius %22.15e \n",radius);
	}
	else if (step < tl->nu3) {
	  // Reasonable step was taken; leave radius alone
	  if (tl->omega3 < 1.0) {
	    radius = tl->omega3 * TaoMin(norm_d, radius);
	  }
	  else if (tl->omega3 > 1.0) {
	    radius = TaoMax(tl->omega3 * norm_d, radius);
	  }
          info=PetscInfo1(tao,"Reasonable step was taken; leave radius alone %22.15e \n",radius);
	}
	else if (step < tl->nu4) {
	  // Full step taken; increase the radius
	  radius = TaoMax(tl->omega4 * norm_d, radius);
          info=PetscInfo1(tao,"Full step taken; increase the radius %22.15e \n",radius);
	}
	else {
	  // More than full step taken; increase the radius
	  radius = TaoMax(tl->omega5 * norm_d, radius);
          info=PetscInfo1(tao,"More than full step taken; increase the radius %22.15e \n",radius);
	}
      }
      else {
	// Newton step was not good; reduce the radius
	radius = tl->omega1 * TaoMin(norm_d, radius);
        info=PetscInfo1(tao,"Newton step was not good; reduce the radius %22.15e \n",radius);
      }
    }
    else {
      // Trust-region step is accepted
      info=PetscInfo(tao,"trust-region step accepted \n");
      info = X->CopyFrom(W); CHKERRQ(info);
      f = ftrial;
      info = TaoComputeGradient(tao, X, G); CHKERRQ(info);
      ++tl->trust;
    }

    // The radius may have been increased; modify if it is too large
    radius = TaoMin(radius, tl->max_radius);

    // Check for termination
    info = G->Norm2(&gnorm); CHKERRQ(info);
    if (TaoInfOrNaN(f) || TaoInfOrNaN(gnorm)) {
      SETERRQ(1,"User provided compute function generated Not-a-Number");
    }
    needH = 1;
    
    info = TaoMonitor(tao, iter, f, gnorm, 0.0, radius, &reason); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_NTL"
static int TaoSetUp_NTL(TAO_SOLVER tao, void *solver)
{
  TAO_NTL *tl = (TAO_NTL *)solver;
  TaoVec *X;
  TaoMat *H;
  int info;

  TaoFunctionBegin;

  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = X->Clone(&tl->G); CHKERRQ(info);
  info = X->Clone(&tl->D); CHKERRQ(info);
  info = X->Clone(&tl->W); CHKERRQ(info);
  info = X->Clone(&tl->XL); CHKERRQ(info);
  info = X->Clone(&tl->XU); CHKERRQ(info);
  info = X->Clone(&tl->PG); CHKERRQ(info);

  info = X->Clone(&tl->Xold); CHKERRQ(info);
  info = X->Clone(&tl->Gold); CHKERRQ(info);

  tl->Diag = 0;
  tl->M = 0;

  info = TaoSetLagrangianGradientVector(tao, tl->G); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, tl->D); CHKERRQ(info);
  info = TaoSetVariableBounds(tao,tl->XL,tl->XU);CHKERRQ(info);

  // Set linear solver to default for symmetric matrices
  info = TaoGetHessian(tao, &H); CHKERRQ(info);
  info = TaoCreateLinearSolver(tao, H, 200, 0); CHKERRQ(info);

  // Check sizes for compatability
  info = TaoCheckFGH(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_NTL"
static int TaoDestroy_NTL(TAO_SOLVER tao, void *solver)
{
  TAO_NTL *tl = (TAO_NTL *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoVecDestroy(tl->G); CHKERRQ(info);
  info = TaoVecDestroy(tl->D); CHKERRQ(info);
  info = TaoVecDestroy(tl->W); CHKERRQ(info);
  info = TaoVecDestroy(tl->XL); CHKERRQ(info);
  info = TaoVecDestroy(tl->XU); CHKERRQ(info);
  info = TaoVecDestroy(tl->PG); CHKERRQ(info);

  info = TaoVecDestroy(tl->Xold); CHKERRQ(info);
  info = TaoVecDestroy(tl->Gold); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao, 0); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, 0); CHKERRQ(info);

  if (tl->Diag) {
    info = TaoVecDestroy(tl->Diag); CHKERRQ(info);
    tl->Diag = 0;
  }

  if (tl->M) {
    info = TaoMatDestroy(tl->M); CHKERRQ(info);
    tl->M = 0;
  }
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_NTL"
static int TaoSetOptions_NTL(TAO_SOLVER tao, void *solver)
{
  TAO_NTL *tl = (TAO_NTL *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoOptionsHead("Newton line search method for unconstrained optimization"); CHKERRQ(info);
  info = TaoOptionList("-tao_ntl_ksp_type", "ksp type", "", NTL_KSP, NTL_KSP_TYPES, NTL_KSP[tl->ksp_type], &tl->ksp_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_ntl_pc_type", "pc type", "", NTL_PC, NTL_PC_TYPES, NTL_PC[tl->pc_type], &tl->pc_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_ntl_bfgs_scale_type", "bfgs scale type", "", BFGS_SCALE, BFGS_SCALE_TYPES, BFGS_SCALE[tl->bfgs_scale_type], &tl->bfgs_scale_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_ntl_init_type", "radius initialization type", "", NTL_INIT, NTL_INIT_TYPES, NTL_INIT[tl->init_type], &tl->init_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_ntl_update_type", "radius update type", "", NTL_UPDATE, NTL_UPDATE_TYPES, NTL_UPDATE[tl->update_type], &tl->update_type, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_eta1", "poor steplength; reduce radius", "", tl->eta1, &tl->eta1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_eta2", "reasonable steplength; leave radius alone", "", tl->eta2, &tl->eta2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_eta3", "good steplength; increase radius", "", tl->eta3, &tl->eta3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_eta4", "excellent steplength; greatly increase radius", "", tl->eta4, &tl->eta4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_alpha1", "", "", tl->alpha1, &tl->alpha1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_alpha2", "", "", tl->alpha2, &tl->alpha2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_alpha3", "", "", tl->alpha3, &tl->alpha3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_alpha4", "", "", tl->alpha4, &tl->alpha4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_alpha5", "", "", tl->alpha5, &tl->alpha5, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_nu1", "poor steplength; reduce radius", "", tl->nu1, &tl->nu1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_nu2", "reasonable steplength; leave radius alone", "", tl->nu2, &tl->nu2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_nu3", "good steplength; increase radius", "", tl->nu3, &tl->nu3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_nu4", "excellent steplength; greatly increase radius", "", tl->nu4, &tl->nu4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_omega1", "", "", tl->omega1, &tl->omega1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_omega2", "", "", tl->omega2, &tl->omega2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_omega3", "", "", tl->omega3, &tl->omega3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_omega4", "", "", tl->omega4, &tl->omega4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_omega5", "", "", tl->omega5, &tl->omega5, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_mu1_i", "", "", tl->mu1_i, &tl->mu1_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_mu2_i", "", "", tl->mu2_i, &tl->mu2_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_gamma1_i", "", "", tl->gamma1_i, &tl->gamma1_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_gamma2_i", "", "", tl->gamma2_i, &tl->gamma2_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_gamma3_i", "", "", tl->gamma3_i, &tl->gamma3_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_gamma4_i", "", "", tl->gamma4_i, &tl->gamma4_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_theta_i", "", "", tl->theta_i, &tl->theta_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_mu1", "", "", tl->mu1, &tl->mu1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_mu2", "", "", tl->mu2, &tl->mu2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_gamma1", "", "", tl->gamma1, &tl->gamma1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_gamma2", "", "", tl->gamma2, &tl->gamma2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_gamma3", "", "", tl->gamma3, &tl->gamma3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_gamma4", "", "", tl->gamma4, &tl->gamma4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_theta", "", "", tl->theta, &tl->theta, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_min_radius", "lower bound on initial radius", "", tl->min_radius, &tl->min_radius, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_max_radius", "upper bound on radius", "", tl->max_radius, &tl->max_radius, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_ntl_epsilon", "tolerance used when computing actual and predicted reduction", "", tl->epsilon, &tl->epsilon, 0); CHKERRQ(info);
  info = TaoLineSearchSetFromOptions(tao); CHKERRQ(info);
  info = TaoOptionsTail(); CHKERRQ(info);
  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_NTL"
static int TaoView_NTL(TAO_SOLVER tao,void* solver)
{
  TAO_NTL *tl = (TAO_NTL *)solver;
  int info;

  TaoFunctionBegin;
  if (NTL_PC_BFGS == tl->pc_type && tl->M) {
    info = TaoPrintInt(tao, "  Rejected matrix updates: %d\n", tl->M->GetRejects()); CHKERRQ(info);
  }
  info = TaoPrintInt(tao, "  Trust-region steps: %d\n", tl->trust); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Newton search steps: %d\n", tl->newt); CHKERRQ(info);
  info = TaoPrintInt(tao, "  BFGS search steps: %d\n", tl->bfgs); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Scaled gradient search steps: %d\n", tl->sgrad); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Gradient search steps: %d\n", tl->grad); CHKERRQ(info);
  info = TaoLineSearchView(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_NTL"
int TaoCreate_NTL(TAO_SOLVER tao)
{
  TAO_NTL *tl;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_NTL, &tl); CHKERRQ(info);
  info = PetscLogObjectMemory(tao, sizeof(TAO_NTL)); CHKERRQ(info);

  info = TaoSetTaoSolveRoutine(tao, TaoSolve_NTL, (void *)tl); CHKERRQ(info);
  info = TaoSetTaoSetUpDownRoutines(tao, TaoSetUp_NTL, TaoDestroy_NTL); CHKERRQ(info);
  info = TaoSetTaoOptionsRoutine(tao, TaoSetOptions_NTL); CHKERRQ(info);
  info = TaoSetTaoViewRoutine(tao, TaoView_NTL); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao, 50); CHKERRQ(info);
  info = TaoSetTolerances(tao, 1e-10, 1e-10, 0, 0); CHKERRQ(info);

  info = TaoSetTrustRegionRadius(tao, 100.0); CHKERRQ(info);
  info = TaoSetTrustRegionTolerance(tao, 1.0e-12); CHKERRQ(info);

  // Default values for trust-region radius update based on steplength
  tl->nu1 = 0.25;
  tl->nu2 = 0.50;
  tl->nu3 = 1.00;
  tl->nu4 = 1.25;

  tl->omega1 = 0.25;
  tl->omega2 = 0.50;
  tl->omega3 = 1.00;
  tl->omega4 = 2.00;
  tl->omega5 = 4.00;

  // Default values for trust-region radius update based on reduction
  tl->eta1 = 1.0e-4;
  tl->eta2 = 0.25;
  tl->eta3 = 0.50;
  tl->eta4 = 0.90;

  tl->alpha1 = 0.25;
  tl->alpha2 = 0.50;
  tl->alpha3 = 1.00;
  tl->alpha4 = 2.00;
  tl->alpha5 = 4.00;

  // Default values for trust-region radius update based on interpolation
  tl->mu1 = 0.10;
  tl->mu2 = 0.50;

  tl->gamma1 = 0.25;
  tl->gamma2 = 0.50;
  tl->gamma3 = 2.00;
  tl->gamma4 = 4.00;

  tl->theta = 0.05;

  // Default values for trust region initialization based on interpolation
  tl->mu1_i = 0.35;
  tl->mu2_i = 0.50;

  tl->gamma1_i = 0.0625;
  tl->gamma2_i = 0.5;
  tl->gamma3_i = 2.0;
  tl->gamma4_i = 5.0;
  
  tl->theta_i = 0.25;

  // Remaining parameters
  tl->min_radius = 1.0e-10;
  tl->max_radius = 1.0e10;
  tl->epsilon = 1.0e-6;

  tl->ksp_type        = NTL_KSP_STCG;
  tl->pc_type         = NTL_PC_BFGS;
  tl->bfgs_scale_type = BFGS_SCALE_AHESS;
  tl->init_type       = NTL_INIT_INTERPOLATION;
  tl->update_type     = NTL_UPDATE_REDUCTION;

  info = TaoCreateMoreThuenteBoundLineSearch(tao,0,0.9); CHKERRQ(info);
  TaoFunctionReturn(0);
}
EXTERN_C_END

#endif

