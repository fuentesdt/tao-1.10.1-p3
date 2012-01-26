/*$Id$*/

#include "nls.h"             /*I "tao_solver.h" I*/

#ifdef TAO_USE_PETSC
#include "petscksp.h"
#include "petscpc.h"
#include "src/petsctao/linearsolver/taolinearsolver_petsc.h"
#include "src/petsctao/vector/taovec_petsc.h"

#include "private/kspimpl.h"
#include "private/pcimpl.h"

#define NLS_KSP_CG	0
#define NLS_KSP_NASH	1
#define NLS_KSP_STCG	2
#define NLS_KSP_GLTR	3
#define NLS_KSP_PETSC	4
#define NLS_KSP_TYPES	5

#define NLS_PC_NONE	0
#define NLS_PC_AHESS	1
#define NLS_PC_BFGS	2
#define NLS_PC_PETSC	3
#define NLS_PC_TYPES	4

#define BFGS_SCALE_AHESS	0
#define BFGS_SCALE_PHESS	1
#define BFGS_SCALE_BFGS		2
#define BFGS_SCALE_TYPES	3

#define NLS_INIT_CONSTANT         0
#define NLS_INIT_DIRECTION        1
#define NLS_INIT_INTERPOLATION    2
#define NLS_INIT_TYPES            3

#define NLS_UPDATE_STEP           0
#define NLS_UPDATE_REDUCTION      1
#define NLS_UPDATE_INTERPOLATION  2
#define NLS_UPDATE_TYPES          3

static const char *NLS_KSP[64] = {
  "cg", "nash", "stcg", "gltr", "petsc"
};

static const char *NLS_PC[64] = {
  "none", "ahess", "bfgs", "petsc"
};

static const char *BFGS_SCALE[64] = {
  "ahess", "phess", "bfgs"
};

static const char *NLS_INIT[64] = {
  "constant", "direction", "interpolation"
};

static const char *NLS_UPDATE[64] = {
  "step", "reduction", "interpolation"
};

// Routine for BFGS preconditioner

#undef __FUNCT__
#define __FUNCT__ "bfgs_apply"
static PetscErrorCode bfgs_apply(PC pc, Vec xin, Vec xout)
{
  TaoLMVMMat *M ;
  TaoVecPetsc Xin(xin);
  TaoVecPetsc Xout(xout);
  TaoTruth info2;
  int info;

  PetscFunctionBegin;

  PetscTruth VerbosePrint = PETSC_FALSE; 
  PetscOptionsGetTruth(PETSC_NULL,"-verboseapp",&VerbosePrint,PETSC_NULL);

  info = PCShellGetContext(pc,(void**)&M); CHKERRQ(info);

  PetscScalar solnNorm,solnDot;
  info = VecNorm(xin,NORM_2,&solnNorm); CHKERRQ(info)
  info=PetscPrintf(PETSC_COMM_WORLD,"bfgs_apply: ||Xin||_2 = %22.15e\n",solnNorm);
  if(VerbosePrint) VecView(xin,0);

  info = M->Solve(&Xin, &Xout, &info2); CHKERRQ(info);

  info = VecNorm(xout,NORM_2,&solnNorm); CHKERRQ(info)
  info = VecDot(xin,xout,&solnDot); CHKERRQ(info)
  info=PetscPrintf(PETSC_COMM_WORLD,"bfgs_apply: ||Xout||_2 = %22.15e, Xin^T Xout= %22.15e\n",solnNorm,solnDot);
  if(VerbosePrint) VecView(xout,0);
  PetscFunctionReturn(0);
}

// Implements Newton's Method with a line search approach for solving
// unconstrained minimization problems.  A More'-Thuente line search 
// is used to guarantee that the bfgs preconditioner remains positive
// definite.
//
// The method can shift the Hessian matrix.  The shifting procedure is
// adapted from the PATH algorithm for solving complementarity
// problems.
//
// The linear system solve should be done with a conjugate gradient
// method, although any method can be used.

#define NLS_NEWTON 		0
#define NLS_BFGS 		1
#define NLS_SCALED_GRADIENT 	2
#define NLS_GRADIENT 		3

#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_NLS"
static int TaoSolve_NLS(TAO_SOLVER tao, void *solver)
{
  TAO_NLS *ls = (TAO_NLS *)solver;
  TaoVec *X, *G = ls->G, *D = ls->D, *W = ls->W, *PG=ls->PG;
  TaoVec *Xold = ls->Xold, *Gold = ls->Gold, *Diag = ls->Diag;
  TaoMat *H;
  TaoLMVMMat *M = ls->M;

  TaoLinearSolver *tls;
  TaoLinearSolverPetsc *pls;

  KSP pksp;
  PC ppc;

  KSPConvergedReason ksp_reason;
  TaoTerminateReason reason;
  TaoTruth success;
  
  double fmin, ftrial, f_full, prered, actred, kappa, sigma;
  double tau, tau_1, tau_2, tau_max, tau_min, max_radius;
  double f, fold, gdx, gnorm, pert;
  double step = 1.0;

  double delta;
  double radius, norm_d = 0.0, e_min;

  int info;
  TaoInt stepType;
  TaoInt iter = 0, status = 0;
  TaoInt bfgsUpdates = 0;
  TaoInt needH;

  TaoInt i_max = 5;
  TaoInt j_max = 1;
  TaoInt i, j;

  TaoFunctionBegin;
  // Initialized variables
  pert = ls->sval;

  ls->ksp_atol = 0;
  ls->ksp_rtol = 0;
  ls->ksp_dtol = 0;
  ls->ksp_ctol = 0;
  ls->ksp_negc = 0;
  ls->ksp_iter = 0;
  ls->ksp_othr = 0;

  // Initialize trust-region radius when using nash, stcg, or gltr
  // Will be reset during the first iteration
  if (NLS_KSP_NASH == ls->ksp_type ||
      NLS_KSP_STCG == ls->ksp_type || 
      NLS_KSP_GLTR == ls->ksp_type) {
    info = TaoGetInitialTrustRegionRadius(tao, &radius); CHKERRQ(info);
    if (radius < 0.0) {
      SETERRQ(1, "Initial radius negative");
    }

    // Modify the radius if it is too large or small
    radius = TaoMax(radius, ls->min_radius);
    radius = TaoMin(radius, ls->max_radius);
  }

  // Get vectors we will need
  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = TaoGetHessian(tao, &H); CHKERRQ(info);

  /*   Project the current point onto the feasible set */
  TaoVec *XU, *XL;
  info = TaoGetVariableBounds(tao,&XL,&XU);CHKERRQ(info);
  info = TaoEvaluateVariableBounds(tao,XL,XU); CHKERRQ(info);
  info = X->Median(XL,X,XU); CHKERRQ(info);
  
  if (NLS_PC_BFGS == ls->pc_type && !M) {
    ls->M = new TaoLMVMMat(X);
    M = ls->M;
  }

  // possibly use BFGS as initial guess
  PetscTruth GradientInitialGuess = PETSC_FALSE;
  PetscOptionsGetTruth(PETSC_NULL,"-tao_nls_gradient_initial_guess",&GradientInitialGuess,PETSC_NULL);
  TaoLMVMMat *InitGuessLMVM  = 0;
  if(!GradientInitialGuess) InitGuessLMVM = new TaoLMVMMat(X);

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
  if ((NLS_PC_BFGS == ls->pc_type) && 
      (BFGS_SCALE_BFGS != ls->bfgs_scale_type)) {
    if (!Diag) {
      info = X->Clone(&ls->Diag); CHKERRQ(info);
      Diag = ls->Diag;
    }
  }

  // Modify the linear solver to a conjugate gradient method
  info = TaoGetLinearSolver(tao, &tls); CHKERRQ(info);
  pls  = dynamic_cast <TaoLinearSolverPetsc *> (tls);

  pksp = pls->GetKSP();
  switch(ls->ksp_type) {
  case NLS_KSP_CG:
    info = KSPSetType(pksp, KSPCG); CHKERRQ(info);
    if (pksp->ops->setfromoptions) {
      (*pksp->ops->setfromoptions)(pksp);
    }
    break;

  case NLS_KSP_NASH:
    info = KSPSetType(pksp, KSPNASH); CHKERRQ(info);
    if (pksp->ops->setfromoptions) {
      (*pksp->ops->setfromoptions)(pksp);
    }
    break;

  case NLS_KSP_STCG:
    info = KSPSetType(pksp, KSPSTCG); CHKERRQ(info);
    if (pksp->ops->setfromoptions) {
      (*pksp->ops->setfromoptions)(pksp);
    }
    break;

  case NLS_KSP_GLTR:
    info = KSPSetType(pksp, KSPGLTR); CHKERRQ(info);
    if (pksp->ops->setfromoptions) {
      (*pksp->ops->setfromoptions)(pksp);
    }
    break;

  default:
    // Use the method set by the ksp_type
    break;
  }

  // Modify the preconditioner to use the bfgs approximation
  info = KSPGetPC(pksp, &ppc); CHKERRQ(info);
  switch(ls->pc_type) {
  case NLS_PC_NONE:
    info = PCSetType(ppc, PCNONE); CHKERRQ(info);
    if (ppc->ops->setfromoptions) {
      (*ppc->ops->setfromoptions)(ppc);
    }
    break;

  case NLS_PC_AHESS:
    info = PCSetType(ppc, PCJACOBI); CHKERRQ(info);
    if (ppc->ops->setfromoptions) {
      (*ppc->ops->setfromoptions)(ppc);
    }
    info = PCJacobiSetUseAbs(ppc); CHKERRQ(info);
    break;

  case NLS_PC_BFGS:
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
  if (NLS_KSP_NASH == ls->ksp_type ||
      NLS_KSP_STCG == ls->ksp_type || 
      NLS_KSP_GLTR == ls->ksp_type) {
    switch(ls->init_type) {
    case NLS_INIT_CONSTANT:
      // Use the initial radius specified
      break;

    case NLS_INIT_INTERPOLATION:
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
            tau = ls->gamma1_i;
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
            if ((fabs(actred) <= ls->epsilon) && 
                (fabs(prered) <= ls->epsilon)) {
              kappa = 1.0;
            }
            else {
              kappa = actred / prered;
            }
  
            tau_1 = ls->theta_i * gnorm * radius / (ls->theta_i * gnorm * radius + (1.0 - ls->theta_i) * prered - actred);
            tau_2 = ls->theta_i * gnorm * radius / (ls->theta_i * gnorm * radius - (1.0 + ls->theta_i) * prered + actred);
            tau_min = TaoMin(tau_1, tau_2);
            tau_max = TaoMax(tau_1, tau_2);
  
            if (fabs(kappa - 1.0) <= ls->mu1_i) {
              // Great agreement
              max_radius = TaoMax(max_radius, radius);
  
              if (tau_max < 1.0) {
                tau = ls->gamma3_i;
              }
              else if (tau_max > ls->gamma4_i) {
                tau = ls->gamma4_i;
              }
              else if (tau_1 >= 1.0 && tau_1 <= ls->gamma4_i && tau_2 < 1.0) {
                tau = tau_1;
              }
              else if (tau_2 >= 1.0 && tau_2 <= ls->gamma4_i && tau_1 < 1.0) {
                tau = tau_2;
              }
              else {
                tau = tau_max;
              }
            }
            else if (fabs(kappa - 1.0) <= ls->mu2_i) {
              // Good agreement
              max_radius = TaoMax(max_radius, radius);
  
              if (tau_max < ls->gamma2_i) {
                tau = ls->gamma2_i;
              }
              else if (tau_max > ls->gamma3_i) {
                tau = ls->gamma3_i;
              }
              else {
                tau = tau_max;
              }
            }
            else {
              // Not good agreement
              if (tau_min > 1.0) {
                tau = ls->gamma2_i;
              }
              else if (tau_max < ls->gamma1_i) {
                tau = ls->gamma1_i;
              }
              else if ((tau_min < ls->gamma1_i) && (tau_max >= 1.0)) {
                tau = ls->gamma1_i;
              }
              else if ((tau_1 >= ls->gamma1_i) && (tau_1 < 1.0) &&
                       ((tau_2 < ls->gamma1_i) || (tau_2 >= 1.0))) {
                tau = tau_1;
              }
              else if ((tau_2 >= ls->gamma1_i) && (tau_2 < 1.0) &&
                       ((tau_1 < ls->gamma1_i) || (tau_2 >= 1.0))) {
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
      radius = TaoMax(radius, ls->min_radius);
      radius = TaoMin(radius, ls->max_radius);
      break;

    default:
      // Norm of the first direction will initialize radius
      radius = 0.0;
      break;
    }
  } 

  // Set initial scaling for the BFGS preconditioner
  // This step is done after computing the initial trust-region radius
  // since the function value may have decreased
  if (NLS_PC_BFGS == ls->pc_type) {
    if (f != 0.0) {
      delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
    }
    else {
      delta = 2.0 / (gnorm*gnorm);
    }
    info = M->SetDelta(delta); CHKERRQ(info);
    info = M->Reset(); CHKERRQ(info);
  }

  if(InitGuessLMVM) {
    if (f != 0.0) {
      delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
    }
    else {
      delta = 2.0 / (gnorm*gnorm);
    }
    info = InitGuessLMVM->SetDelta(delta); CHKERRQ(info);
    info = InitGuessLMVM->Reset(); CHKERRQ(info);
  }

  // Set counter for gradient/reset steps
  ls->newt = 0;
  ls->bfgs = 0;
  ls->sgrad = 0;
  ls->grad = 0;

  // Have not converged; continue with Newton method
  while (reason == TAO_CONTINUE_ITERATING) {
    ++iter;

    // Compute the Hessian
    if (needH) {
      info = TaoComputeHessian(tao, X, H); CHKERRQ(info);
      needH = 0;
    }

    if ((NLS_PC_BFGS == ls->pc_type) && 
        (BFGS_SCALE_AHESS == ls->bfgs_scale_type)) {
      // Obtain diagonal for the bfgs preconditioner 
      info = H->GetDiagonal(Diag); CHKERRQ(info);
      info = Diag->AbsoluteValue(); CHKERRQ(info);
      info = Diag->Reciprocal(); CHKERRQ(info);
      info = M->SetScale(Diag); CHKERRQ(info);
      M->View();
    }
 
    // Shift the Hessian matrix
    info=PetscInfo1(tao,"TaoSolve_NLS: pert %22.15e \n",pert);
    if (pert > 0) {
      info = H->ShiftDiagonal(pert); CHKERRQ(info);
    }
    
    if (NLS_PC_BFGS == ls->pc_type) {
      if (BFGS_SCALE_PHESS == ls->bfgs_scale_type) {
        // Obtain diagonal for the bfgs preconditioner 
        info = H->GetDiagonal(Diag); CHKERRQ(info);
        info = Diag->AbsoluteValue(); CHKERRQ(info);
        info = Diag->Reciprocal(); CHKERRQ(info);
        info = M->SetScale(Diag); CHKERRQ(info);
      }

      // Update the limited memory preconditioner
      info = M->Update(X, G); CHKERRQ(info);
      M->View();
      ++bfgsUpdates;
    }

    if(InitGuessLMVM) {
      // Update the limited memory initial guess
      info = InitGuessLMVM->Update(X, G); CHKERRQ(info);
      InitGuessLMVM->View();
    }

    // use constrained norm for tolerance
    PetscScalar boundNorm;
    info = PG->BoundGradientProjection(G,XL,X,XU);CHKERRQ(info);
    info = PG->Norm2(&boundNorm); CHKERRQ(info);
    
    PetscScalar ewAtol  = PetscMin(0.5,boundNorm)*boundNorm;
    info = KSPSetTolerances(pksp,PETSC_DEFAULT,ewAtol,
                            PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(info);
    info=PetscInfo2(tao,"TaoSolve_NLS: gnorm=%22.12e, boundNorm=%22.12e\n",gnorm,boundNorm);
    pksp->printreason = PETSC_TRUE;
    info = KSPView(pksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(info);

    // use gradient as initial guess
    // use bfgs as initial guess
    if( InitGuessLMVM ) 
     {
      info=PetscInfo(tao,"TaoSolve_NLS: use bfgs init guess \n");
      info = InitGuessLMVM->Solve(G, D, &success);
     }
    else
      info = D->CopyFrom(G); 
    CHKERRQ(info);

    // Solve the Newton system of equations
    info = TaoPreLinearSolve(tao, H); CHKERRQ(info);
    if (NLS_KSP_NASH == ls->ksp_type ||
        NLS_KSP_STCG == ls->ksp_type || 
        NLS_KSP_GLTR == ls->ksp_type) {
      info = TaoLinearSolveTrustRegion(tao, H, G, D, radius, &success); CHKERRQ(info);
      info = pls->GetNormDirection(&norm_d); CHKERRQ(info);
      if (0.0 == radius) {
        // Radius was uninitialized; use the norm of the direction
        if (norm_d > 0.0) {
          radius = norm_d;

          // Modify the radius if it is too large or small
          radius = TaoMax(radius, ls->min_radius);
          radius = TaoMin(radius, ls->max_radius);
        }
        else {
          // The direction was bad; set radius to default value and re-solve 
	  // the trust-region subproblem to get a direction
          info = TaoGetInitialTrustRegionRadius(tao, &radius); CHKERRQ(info);

          // Modify the radius if it is too large or small
          radius = TaoMax(radius, ls->min_radius);
          radius = TaoMin(radius, ls->max_radius);

          info = TaoLinearSolveTrustRegion(tao, H, G, D, radius, &success); CHKERRQ(info);
          info = pls->GetNormDirection(&norm_d); CHKERRQ(info);
          if (norm_d == 0.0) {
            SETERRQ(1, "Initial direction zero");
          }
        }
      }
    }
    else {
      info = TaoLinearSolve(tao, H, G, D, &success); CHKERRQ(info);
    }
    info = D->Negate(); CHKERRQ(info);

    info = KSPGetConvergedReason(pksp, &ksp_reason); CHKERRQ(info);
    if ((KSP_DIVERGED_INDEFINITE_PC == ksp_reason) &&
        (NLS_PC_BFGS == ls->pc_type) && (bfgsUpdates > 1)) {
      // Preconditioner is numerically indefinite; reset the
      // approximate if using BFGS preconditioning.

      info=PetscInfo(tao,"TaoSolve_NLS: Indefinite_PC reset BFGS\n");
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

    if (KSP_CONVERGED_ATOL == ksp_reason) {
      ++ls->ksp_atol;
    }
    else if (KSP_CONVERGED_RTOL == ksp_reason) {
      ++ls->ksp_rtol;
    }
    else if (KSP_CONVERGED_CG_CONSTRAINED == ksp_reason) {
      ++ls->ksp_ctol;
    }
    else if (KSP_CONVERGED_CG_NEG_CURVE == ksp_reason) {
      ++ls->ksp_negc;
    }
    else if (KSP_DIVERGED_DTOL == ksp_reason) {
      ++ls->ksp_dtol;
    }
    else if (KSP_DIVERGED_ITS == ksp_reason) {
      ++ls->ksp_iter;
    }
    else {
      ++ls->ksp_othr;
    } 

    // Check for success (descent direction)
    info = D->Dot(G, &gdx); CHKERRQ(info);
    if ((gdx >= 0.0) || TaoInfOrNaN(gdx)) {
      // Newton step is not descent or direction produced Inf or NaN
      // Update the perturbation for next time
      info=PetscInfo2(tao,"TaoSolve_NLS: Newton step is not descent or direction produced Inf or NaN,pert %22.12e, gdx %22.12e \n", pert,gdx);
      if (pert <= 0.0) {
	// Initialize the perturbation
	pert = TaoMin(ls->imax, TaoMax(ls->imin, ls->imfac * gnorm));
        if (NLS_KSP_GLTR == ls->ksp_type) {
          info = pls->GetMinEig(&e_min); CHKERRQ(info);
	  pert = TaoMax(pert, -e_min);
        }
      }
      else {
	// Increase the perturbation
	pert = TaoMin(ls->pmax, TaoMax(ls->pgfac * pert, ls->pmgfac * gnorm));
      }

      if(InitGuessLMVM) {
        info=PetscInfo(tao,"TaoSolve_NLS: reset init guess \n");
        if (f != 0.0) {
          delta = 2.0 * TaoAbsDouble(f) / (gnorm*gnorm);
        }
        else {
          delta = 2.0 / (gnorm*gnorm);
        }
        info = InitGuessLMVM->SetDelta(delta); CHKERRQ(info);
        info = InitGuessLMVM->Reset(); CHKERRQ(info);
      }

      if (NLS_PC_BFGS != ls->pc_type) {
	// We don't have the bfgs matrix around and updated
        // Must use gradient direction in this case
        info=PetscInfo(tao,"TaoSolve_NLS: using steepest descent \n");
        info = D->CopyFrom(G); CHKERRQ(info);
        info = D->Negate(); CHKERRQ(info);

	++ls->grad;
        stepType = NLS_GRADIENT;
      }
      else {
        // Attempt to use the BFGS direction
        info=PetscInfo(tao,"TaoSolve_NLS: trying BFGS Direction\n");
        info = M->Solve(G, D, &success); CHKERRQ(info);
        info = D->Negate(); CHKERRQ(info);

        // Check for success (descent direction)
        info = D->Dot(G, &gdx); CHKERRQ(info);
        if ((gdx >= 0) || TaoInfOrNaN(gdx)) {
          // BFGS direction is not descent or direction produced not a number
          // We can assert bfgsUpdates > 1 in this case because
          // the first solve produces the scaled gradient direction,
          // which is guaranteed to be descent
	  //
          // Use steepest descent direction (scaled)
          info=PetscInfo(tao,"TaoSolve_NLS: BFGS Direction fail\n");

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
          ++ls->sgrad;
          stepType = NLS_SCALED_GRADIENT;
        }
        else {
          if (1 == bfgsUpdates) {
	    // The first BFGS direction is always the scaled gradient
            ++ls->sgrad;
            stepType = NLS_SCALED_GRADIENT;
          }
          else {
            ++ls->bfgs;
            stepType = NLS_BFGS;
          }
        }
      }
    }
    else {
      // Computed Newton step is descent
      info=PetscInfo(tao,"TaoSolve_NLS: Newton step is descent\n");
      switch (ksp_reason) {
      case KSP_DIVERGED_NAN:
      case KSP_DIVERGED_BREAKDOWN:
      case KSP_DIVERGED_INDEFINITE_MAT:
      case KSP_DIVERGED_INDEFINITE_PC:
      case KSP_CONVERGED_CG_NEG_CURVE:
        // Matrix or preconditioner is indefinite; increase perturbation
        if (pert <= 0.0) {
	  // Initialize the perturbation
          pert = TaoMin(ls->imax, TaoMax(ls->imin, ls->imfac * gnorm));
          if (NLS_KSP_GLTR == ls->ksp_type) {
            info = pls->GetMinEig(&e_min); CHKERRQ(info);
	    pert = TaoMax(pert, -e_min);
          }
        }
        else {
	  // Increase the perturbation
	  pert = TaoMin(ls->pmax, TaoMax(ls->pgfac * pert, ls->pmgfac * gnorm));
        }
        break;

      default:
        // Newton step computation is good; decrease perturbation
        pert = TaoMin(ls->psfac * pert, ls->pmsfac * gnorm);
        if (pert < ls->pmin) {
	  pert = 0.0;
        }
        break; 
      }

      ++ls->newt;
      stepType = NLS_NEWTON;
    }

    // Perform the linesearch
    fold = f;
    info = Xold->CopyFrom(X); CHKERRQ(info);
    info = Gold->CopyFrom(G); CHKERRQ(info);

    step = 1.0;
    info = TaoLineSearchApply(tao, X, G, D, W, &f, &f_full, &step, &status); CHKERRQ(info);
    info=PetscInfo1(tao,"TaoSolve_NLS: stepType %d \n",stepType);

    while (status && stepType != NLS_GRADIENT) {
      // Linesearch failed
      f = fold;
      info = X->CopyFrom(Xold); CHKERRQ(info);
      info = G->CopyFrom(Gold); CHKERRQ(info);

      switch(stepType) {
      case NLS_NEWTON:
        // Failed to obtain acceptable iterate with Newton step
        // Update the perturbation for next time
        if (pert <= 0.0) {
          // Initialize the perturbation
          pert = TaoMin(ls->imax, TaoMax(ls->imin, ls->imfac * gnorm));
          if (NLS_KSP_GLTR == ls->ksp_type) {
            info = pls->GetMinEig(&e_min); CHKERRQ(info);
	    pert = TaoMax(pert, -e_min);
          }
        }
        else {
          // Increase the perturbation
          pert = TaoMin(ls->pmax, TaoMax(ls->pgfac * pert, ls->pmgfac * gnorm));
        }

        if (NLS_PC_BFGS != ls->pc_type) {
	  // We don't have the bfgs matrix around and being updated
          // Must use gradient direction in this case
          info = D->CopyFrom(G); CHKERRQ(info);

	  ++ls->grad;
          stepType = NLS_GRADIENT;
        }
        else {
          // Attempt to use the BFGS direction
          info = M->Solve(G, D, &success); CHKERRQ(info);

          // Check for success (descent direction)
          info = D->Dot(G, &gdx); CHKERRQ(info);
          if ((gdx <= 0) || TaoInfOrNaN(gdx)) {
            // BFGS direction is not descent or direction produced not a number
            // We can assert bfgsUpdates > 1 in this case
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
            ++ls->sgrad;
            stepType = NLS_SCALED_GRADIENT;
          }
          else {
            if (1 == bfgsUpdates) {
	      // The first BFGS direction is always the scaled gradient
              ++ls->sgrad;
              stepType = NLS_SCALED_GRADIENT;
            }
            else {
              ++ls->bfgs;
              stepType = NLS_BFGS;
            }
          }
        }
	break;

      case NLS_BFGS:
        // Can only enter if pc_type == NLS_PC_BFGS
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
        ++ls->sgrad;
        stepType = NLS_SCALED_GRADIENT;
        break;

      case NLS_SCALED_GRADIENT:
        // Can only enter if pc_type == NLS_PC_BFGS
        // The scaled gradient step did not produce a new iterate;
        // attemp to use the gradient direction.
        // Need to make sure we are not using a different diagonal scaling
	info = M->SetScale(0); CHKERRQ(info);
        info = M->SetDelta(1.0); CHKERRQ(info);
        info = M->Reset(); CHKERRQ(info);
        info = M->Update(X, G); CHKERRQ(info);
        info = M->Solve(G, D, &success); CHKERRQ(info);

        bfgsUpdates = 1;
	++ls->grad;
        stepType = NLS_GRADIENT;
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
      step = 0.0;
    }

    // Update trust region radius
    if (NLS_KSP_NASH == ls->ksp_type ||
        NLS_KSP_STCG == ls->ksp_type || 
        NLS_KSP_GLTR == ls->ksp_type) {
      switch(ls->update_type) {
      case NLS_UPDATE_STEP:
        if (stepType == NLS_NEWTON) {
          if (step < ls->nu1) {
            // Very bad step taken; reduce radius
            radius = ls->omega1 * TaoMin(norm_d, radius);
          }
          else if (step < ls->nu2) {
            // Reasonably bad step taken; reduce radius
            radius = ls->omega2 * TaoMin(norm_d, radius);
          }
          else if (step < ls->nu3) {
            // Reasonable step was taken; leave radius alone
            if (ls->omega3 < 1.0) {
              radius = ls->omega3 * TaoMin(norm_d, radius);
            }
            else if (ls->omega3 > 1.0) {
              radius = TaoMax(ls->omega3 * norm_d, radius);  
            }
          }
          else if (step < ls->nu4) {
            // Full step taken; increase the radius
            radius = TaoMax(ls->omega4 * norm_d, radius);  
          }
          else {
            // More than full step taken; increase the radius
            radius = TaoMax(ls->omega5 * norm_d, radius);  
          }
        }
        else {
          // Newton step was not good; reduce the radius
          radius = ls->omega1 * TaoMin(norm_d, radius);
        }
        break;

      case NLS_UPDATE_REDUCTION:
        if (stepType == NLS_NEWTON) {
	  // Get predicted reduction
          info = pls->GetObjFcn(&prered); CHKERRQ(info);

          if (prered >= 0.0) {
            // The predicted reduction has the wrong sign.  This cannot
            // happen in infinite precision arithmetic.  Step should
            // be rejected!
            radius = ls->alpha1 * TaoMin(radius, norm_d);
          }
          else {
            if (TaoInfOrNaN(f_full)) {
              radius = ls->alpha1 * TaoMin(radius, norm_d);
            }
            else {
              // Compute and actual reduction
              actred = fold - f_full;
              prered = -prered;
              if ((fabs(actred) <= ls->epsilon) && 
                  (fabs(prered) <= ls->epsilon)) {
                kappa = 1.0;
              }
              else {
                kappa = actred / prered;
              }
  
              // Accept of reject the step and update radius
              if (kappa < ls->eta1) {
                // Very bad step
                radius = ls->alpha1 * TaoMin(radius, norm_d);
              }
              else if (kappa < ls->eta2) {
                // Marginal bad step
                radius = ls->alpha2 * TaoMin(radius, norm_d);
              }
              else if (kappa < ls->eta3) {
                // Reasonable step
                if (ls->alpha3 < 1.0) {
                  radius = ls->alpha3 * TaoMin(norm_d, radius);
                }
                else if (ls->alpha3 > 1.0) {
                  radius = TaoMax(ls->alpha3 * norm_d, radius);  
                }
              }
              else if (kappa < ls->eta4) {
                // Good step
                radius = TaoMax(ls->alpha4 * norm_d, radius);
              }
              else {
                // Very good step
                radius = TaoMax(ls->alpha5 * norm_d, radius);
              }
            }
          }
        }
        else {
          // Newton step was not good; reduce the radius
          radius = ls->alpha1 * TaoMin(norm_d, radius);
        }
        break;

      default:
        if (stepType == NLS_NEWTON) {
          // Get predicted reduction
          info = pls->GetObjFcn(&prered); CHKERRQ(info);

          if (prered >= 0.0) {
            // The predicted reduction has the wrong sign.  This cannot
            // happen in infinite precision arithmetic.  Step should
            // be rejected!
            radius = ls->gamma1 * TaoMin(radius, norm_d);
          }
          else {
            if (TaoInfOrNaN(f_full)) {
              radius = ls->gamma1 * TaoMin(radius, norm_d);
            }
            else {
              actred = fold - f_full;
              prered = -prered;
              if ((fabs(actred) <= ls->epsilon) && 
                  (fabs(prered) <= ls->epsilon)) {
                kappa = 1.0;
              }
              else {
                kappa = actred / prered;
              }

              tau_1 = ls->theta * gdx / (ls->theta * gdx - (1.0 - ls->theta) * prered + actred);
              tau_2 = ls->theta * gdx / (ls->theta * gdx + (1.0 + ls->theta) * prered - actred);
              tau_min = TaoMin(tau_1, tau_2);
              tau_max = TaoMax(tau_1, tau_2);

              if (kappa >= 1.0 - ls->mu1) {
                // Great agreement
                if (tau_max < 1.0) {
                  radius = TaoMax(radius, ls->gamma3 * norm_d);
                }
                else if (tau_max > ls->gamma4) {
                  radius = TaoMax(radius, ls->gamma4 * norm_d);
                }
                else {
                  radius = TaoMax(radius, tau_max * norm_d);
                }
              }
              else if (kappa >= 1.0 - ls->mu2) {
                // Good agreement

                if (tau_max < ls->gamma2) {
                  radius = ls->gamma2 * TaoMin(radius, norm_d);
                }
                else if (tau_max > ls->gamma3) {
                  radius = TaoMax(radius, ls->gamma3 * norm_d);
                }
                else if (tau_max < 1.0) {
                  radius = tau_max * TaoMin(radius, norm_d);
                }
                else {
                  radius = TaoMax(radius, tau_max * norm_d);
                }
              }
              else {
                // Not good agreement
                if (tau_min > 1.0) {
                  radius = ls->gamma2 * TaoMin(radius, norm_d);
                }
                else if (tau_max < ls->gamma1) {
                  radius = ls->gamma1 * TaoMin(radius, norm_d);
                }
                else if ((tau_min < ls->gamma1) && (tau_max >= 1.0)) {
                  radius = ls->gamma1 * TaoMin(radius, norm_d);
                }
                else if ((tau_1 >= ls->gamma1) && (tau_1 < 1.0) &&
                         ((tau_2 < ls->gamma1) || (tau_2 >= 1.0))) {
                  radius = tau_1 * TaoMin(radius, norm_d);
                }
                else if ((tau_2 >= ls->gamma1) && (tau_2 < 1.0) &&
                         ((tau_1 < ls->gamma1) || (tau_2 >= 1.0))) {
                  radius = tau_2 * TaoMin(radius, norm_d);
                }
                else {
                  radius = tau_max * TaoMin(radius, norm_d);
                }
              }
            } 
          }
        }
        else {
          // Newton step was not good; reduce the radius
          radius = ls->gamma1 * TaoMin(norm_d, radius);
        }
        break;
      }

      // The radius may have been increased; modify if it is too large
      radius = TaoMin(radius, ls->max_radius);
    }

    // Check for termination
    info = G->Norm2(&gnorm); CHKERRQ(info);
    if (TaoInfOrNaN(f) || TaoInfOrNaN(gnorm)) {
      SETERRQ(1,"User provided compute function generated Not-a-Number");
    }
    needH = 1;

    info = TaoMonitor(tao, iter, f, gnorm, 0.0, step, &reason); CHKERRQ(info);
  }

  // clean up
  if(InitGuessLMVM){info=TaoMatDestroy(InitGuessLMVM);CHKERRQ(info);}

  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_NLS"
static int TaoSetUp_NLS(TAO_SOLVER tao, void *solver)
{
  TAO_NLS *ls = (TAO_NLS *)solver;
  TaoVec *X;
  TaoMat *H;
  int info;

  TaoFunctionBegin;

  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = X->Clone(&ls->G); CHKERRQ(info);
  info = X->Clone(&ls->D); CHKERRQ(info);
  info = X->Clone(&ls->W); CHKERRQ(info);
  info = X->Clone(&ls->XL); CHKERRQ(info);
  info = X->Clone(&ls->XU); CHKERRQ(info);
  info = X->Clone(&ls->PG); CHKERRQ(info);

  info = X->Clone(&ls->Xold); CHKERRQ(info);
  info = X->Clone(&ls->Gold); CHKERRQ(info);

  ls->Diag = 0;
  ls->M = 0;

  info = TaoSetLagrangianGradientVector(tao, ls->G); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, ls->D); CHKERRQ(info);
  info = TaoSetVariableBounds(tao,ls->XL,ls->XU);CHKERRQ(info);

  // Set linear solver to default for symmetric matrices
  info = TaoGetHessian(tao, &H); CHKERRQ(info);
  info = TaoCreateLinearSolver(tao, H, 200, 0); CHKERRQ(info);

  // Check sizes for compatability
  info = TaoCheckFGH(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_NLS"
static int TaoDestroy_NLS(TAO_SOLVER tao, void *solver)
{
  TAO_NLS *ls = (TAO_NLS *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoVecDestroy(ls->G); CHKERRQ(info);
  info = TaoVecDestroy(ls->D); CHKERRQ(info);
  info = TaoVecDestroy(ls->W); CHKERRQ(info);
  info = TaoVecDestroy(ls->XL); CHKERRQ(info);
  info = TaoVecDestroy(ls->XU); CHKERRQ(info);
  info = TaoVecDestroy(ls->PG); CHKERRQ(info);

  info = TaoVecDestroy(ls->Xold); CHKERRQ(info);
  info = TaoVecDestroy(ls->Gold); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao, 0); CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao, 0); CHKERRQ(info);

  if (ls->Diag) {
    info = TaoVecDestroy(ls->Diag); CHKERRQ(info);
    ls->Diag = 0;
  }

  if (ls->M) {
    info = TaoMatDestroy(ls->M); CHKERRQ(info);
    ls->M = 0;
  }
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_NLS"
static int TaoSetOptions_NLS(TAO_SOLVER tao, void *solver)
{
  TAO_NLS *ls = (TAO_NLS *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoOptionsHead("Newton line search method for unconstrained optimization"); CHKERRQ(info);
  info = TaoOptionList("-tao_nls_ksp_type", "ksp type", "", NLS_KSP, NLS_KSP_TYPES, NLS_KSP[ls->ksp_type], &ls->ksp_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_nls_pc_type", "pc type", "", NLS_PC, NLS_PC_TYPES, NLS_PC[ls->pc_type], &ls->pc_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_nls_bfgs_scale_type", "bfgs scale type", "", BFGS_SCALE, BFGS_SCALE_TYPES, BFGS_SCALE[ls->bfgs_scale_type], &ls->bfgs_scale_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_nls_init_type", "radius initialization type", "", NLS_INIT, NLS_INIT_TYPES, NLS_INIT[ls->init_type], &ls->init_type, 0); CHKERRQ(info);
  info = TaoOptionList("-tao_nls_update_type", "radius update type", "", NLS_UPDATE, NLS_UPDATE_TYPES, NLS_UPDATE[ls->update_type], &ls->update_type, 0); CHKERRQ(info);
 info = TaoOptionDouble("-tao_nls_sval", "perturbation starting value", "", ls->sval, &ls->sval, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_imin", "minimum initial perturbation", "", ls->imin, &ls->imin, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_imax", "maximum initial perturbation", "", ls->imax, &ls->imax, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_imfac", "initial merit factor", "", ls->imfac, &ls->imfac, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_pmin", "minimum perturbation", "", ls->pmin, &ls->pmin, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_pmax", "maximum perturbation", "", ls->pmax, &ls->pmax, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_pgfac", "growth factor", "", ls->pgfac, &ls->pgfac, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_psfac", "shrink factor", "", ls->psfac, &ls->psfac, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_pmgfac", "merit growth factor", "", ls->pmgfac, &ls->pmgfac, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_pmsfac", "merit shrink factor", "", ls->pmsfac, &ls->pmsfac, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_eta1", "poor steplength; reduce radius", "", ls->eta1, &ls->eta1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_eta2", "reasonable steplength; leave radius alone", "", ls->eta2, &ls->eta2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_eta3", "good steplength; increase radius", "", ls->eta3, &ls->eta3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_eta4", "excellent steplength; greatly increase radius", "", ls->eta4, &ls->eta4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_alpha1", "", "", ls->alpha1, &ls->alpha1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_alpha2", "", "", ls->alpha2, &ls->alpha2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_alpha3", "", "", ls->alpha3, &ls->alpha3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_alpha4", "", "", ls->alpha4, &ls->alpha4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_alpha5", "", "", ls->alpha5, &ls->alpha5, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_nu1", "poor steplength; reduce radius", "", ls->nu1, &ls->nu1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_nu2", "reasonable steplength; leave radius alone", "", ls->nu2, &ls->nu2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_nu3", "good steplength; increase radius", "", ls->nu3, &ls->nu3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_nu4", "excellent steplength; greatly increase radius", "", ls->nu4, &ls->nu4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_omega1", "", "", ls->omega1, &ls->omega1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_omega2", "", "", ls->omega2, &ls->omega2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_omega3", "", "", ls->omega3, &ls->omega3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_omega4", "", "", ls->omega4, &ls->omega4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_omega5", "", "", ls->omega5, &ls->omega5, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_mu1_i", "", "", ls->mu1_i, &ls->mu1_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_mu2_i", "", "", ls->mu2_i, &ls->mu2_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_gamma1_i", "", "", ls->gamma1_i, &ls->gamma1_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_gamma2_i", "", "", ls->gamma2_i, &ls->gamma2_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_gamma3_i", "", "", ls->gamma3_i, &ls->gamma3_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_gamma4_i", "", "", ls->gamma4_i, &ls->gamma4_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_theta_i", "", "", ls->theta_i, &ls->theta_i, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_mu1", "", "", ls->mu1, &ls->mu1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_mu2", "", "", ls->mu2, &ls->mu2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_gamma1", "", "", ls->gamma1, &ls->gamma1, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_gamma2", "", "", ls->gamma2, &ls->gamma2, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_gamma3", "", "", ls->gamma3, &ls->gamma3, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_gamma4", "", "", ls->gamma4, &ls->gamma4, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_theta", "", "", ls->theta, &ls->theta, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_min_radius", "lower bound on initial radius", "", ls->min_radius, &ls->min_radius, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_max_radius", "upper bound on radius", "", ls->max_radius, &ls->max_radius, 0); CHKERRQ(info);
  info = TaoOptionDouble("-tao_nls_epsilon", "tolerance used when computing actual and predicted reduction", "", ls->epsilon, &ls->epsilon, 0); CHKERRQ(info);
  info = TaoLineSearchSetFromOptions(tao); CHKERRQ(info);
  info = TaoOptionsTail(); CHKERRQ(info);
  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_NLS"
static int TaoView_NLS(TAO_SOLVER tao,void* solver)
{
  TAO_NLS *ls = (TAO_NLS *)solver;
  int info;

  TaoFunctionBegin;
  if (NLS_PC_BFGS == ls->pc_type && ls->M) {
    info = TaoPrintInt(tao, "  Rejected matrix updates: %d\n", ls->M->GetRejects()); CHKERRQ(info);
  }
  info = TaoPrintInt(tao, "  Newton steps: %d\n", ls->newt); CHKERRQ(info);
  info = TaoPrintInt(tao, "  BFGS steps: %d\n", ls->bfgs); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Scaled gradient steps: %d\n", ls->sgrad); CHKERRQ(info);
  info = TaoPrintInt(tao, "  Gradient steps: %d\n", ls->grad); CHKERRQ(info);

  info = TaoPrintInt(tao, "  nls ksp atol: %d\n", ls->ksp_atol); CHKERRQ(info);
  info = TaoPrintInt(tao, "  nls ksp rtol: %d\n", ls->ksp_rtol); CHKERRQ(info);
  info = TaoPrintInt(tao, "  nls ksp ctol: %d\n", ls->ksp_ctol); CHKERRQ(info);
  info = TaoPrintInt(tao, "  nls ksp negc: %d\n", ls->ksp_negc); CHKERRQ(info);
  info = TaoPrintInt(tao, "  nls ksp dtol: %d\n", ls->ksp_dtol); CHKERRQ(info);
  info = TaoPrintInt(tao, "  nls ksp iter: %d\n", ls->ksp_iter); CHKERRQ(info);
  info = TaoPrintInt(tao, "  nls ksp othr: %d\n", ls->ksp_othr); CHKERRQ(info);
  info = TaoLineSearchView(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_NLS"
int TaoCreate_NLS(TAO_SOLVER tao)
{
  TAO_NLS *ls;
  int info;

  TaoFunctionBegin;

  info = TaoNew(TAO_NLS, &ls); CHKERRQ(info);
  info = PetscLogObjectMemory(tao, sizeof(TAO_NLS)); CHKERRQ(info);

  info = TaoSetTaoSolveRoutine(tao, TaoSolve_NLS, (void *)ls); CHKERRQ(info);
  info = TaoSetTaoSetUpDownRoutines(tao, TaoSetUp_NLS, TaoDestroy_NLS); CHKERRQ(info);
  info = TaoSetTaoOptionsRoutine(tao, TaoSetOptions_NLS); CHKERRQ(info);
  info = TaoSetTaoViewRoutine(tao, TaoView_NLS); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao, 50); CHKERRQ(info);
  info = TaoSetTolerances(tao, 1e-10, 1e-10, 0, 0); CHKERRQ(info);

  info = TaoSetTrustRegionRadius(tao, 100.0); CHKERRQ(info);
  info = TaoSetTrustRegionTolerance(tao, 1.0e-12); CHKERRQ(info);

  ls->sval   = 0.0;
  ls->imin   = 1.0e-4;
  ls->imax   = 1.0e+2;
  ls->imfac  = 1.0e-1;

  ls->pmin   = 1.0e-12;
  ls->pmax   = 1.0e+2;
  ls->pgfac  = 1.0e+1;
  ls->psfac  = 4.0e-1;
  ls->pmgfac = 1.0e-1;
  ls->pmsfac = 1.0e-1;

  // Default values for trust-region radius update based on steplength
  ls->nu1 = 0.25;
  ls->nu2 = 0.50;
  ls->nu3 = 1.00;
  ls->nu4 = 1.25;

  ls->omega1 = 0.25;
  ls->omega2 = 0.50;
  ls->omega3 = 1.00;
  ls->omega4 = 2.00;
  ls->omega5 = 4.00;

  // Default values for trust-region radius update based on reduction
  ls->eta1 = 1.0e-4;
  ls->eta2 = 0.25;
  ls->eta3 = 0.50;
  ls->eta4 = 0.90;

  ls->alpha1 = 0.25;
  ls->alpha2 = 0.50;
  ls->alpha3 = 1.00;
  ls->alpha4 = 2.00;
  ls->alpha5 = 4.00;

  // Default values for trust-region radius update based on interpolation
  ls->mu1 = 0.10;
  ls->mu2 = 0.50;

  ls->gamma1 = 0.25;
  ls->gamma2 = 0.50;
  ls->gamma3 = 2.00;
  ls->gamma4 = 4.00;

  ls->theta = 0.05;

  // Default values for trust region initialization based on interpolation
  ls->mu1_i = 0.35;
  ls->mu2_i = 0.50;

  ls->gamma1_i = 0.0625;
  ls->gamma2_i = 0.5;
  ls->gamma3_i = 2.0;
  ls->gamma4_i = 5.0;
  
  ls->theta_i = 0.25;

  // Remaining parameters
  ls->min_radius = 1.0e-10;
  ls->max_radius = 1.0e10;
  ls->epsilon = 1.0e-6;

  ls->ksp_type        = NLS_KSP_CG;
  //ls->pc_type         = NLS_PC_BFGS;
  ls->pc_type         = NLS_PC_NONE;
  ls->bfgs_scale_type = BFGS_SCALE_PHESS;
  ls->init_type       = NLS_INIT_INTERPOLATION;
  ls->update_type     = NLS_UPDATE_STEP;

  //info = TaoCreateMoreThuenteLineSearch(tao, 0, 0); CHKERRQ(info);
  info = TaoCreateMoreThuenteBoundLineSearch(tao,0,0.9); CHKERRQ(info);
  TaoFunctionReturn(0);
}
EXTERN_C_END

#endif

