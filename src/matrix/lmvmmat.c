#include <iostream>
#include "lmvmmat.h"
#include "tao_general.h"
#include "tao_solver.h"
#include "taovec.h"
#include "petsc.h"

#define Scale_None		0
#define Scale_Scalar		1
#define Scale_Broyden		2
#define Scale_Types             3

#define Rescale_None		0
#define Rescale_Scalar		1
#define Rescale_GL		2
#define Rescale_Types          	3

#define Limit_None		0
#define Limit_Average		1
#define Limit_Relative		2
#define Limit_Absolute		3
#define Limit_Types		4

#define TAO_ZER_SAFEGUARD	1e-8
#define TAO_INF_SAFEGUARD	1e+8

static const char *Scale_Table[64] = {
  "none", "scalar", "broyden"
};

static const char *Rescale_Table[64] = {
  "none", "scalar", "gl"
};

static const char *Limit_Table[64] = {
  "none", "average", "relative", "absolute"
};

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::TaoLMVMMat"
TaoLMVMMat::TaoLMVMMat(TaoVec *tv)
{
  // Set default values
  lm = 5;
  eps = 0.0;

  limitType = Limit_None;
  scaleType = Scale_Broyden;
  rScaleType = Rescale_Scalar;

  s_alpha = 1.0;
  r_alpha = 1.0;
  r_beta = 0.5;
  mu = 1.0;
  nu = 100.0;

  phi = 0.125;		

  scalar_history = 1;
  rescale_history = 1;

  delta_min = 1e-7;
  delta_max = 100;

  // Begin configuration
  TaoOptionsHead("Limited memory matrix options");
  TaoOptionInt("-tao_lmm_vectors", "vectors to use for approximation", "", lm, &lm, 0);
  TaoOptionDouble("-tao_lmm_limit_mu", "mu limiting factor", "", mu, &mu, 0);
  TaoOptionDouble("-tao_lmm_limit_nu", "nu limiting factor", "", nu, &nu, 0);
  TaoOptionDouble("-tao_lmm_broyden_phi", "phi factor for Broyden scaling", "", phi, &phi, 0);
  TaoOptionDouble("-tao_lmm_scalar_alpha", "alpha factor for scalar scaling", "", s_alpha, &s_alpha, 0);
  TaoOptionDouble("-tao_lmm_rescale_alpha", "alpha factor for rescaling diagonal", "", r_alpha, &r_alpha, 0);
  TaoOptionDouble("-tao_lmm_rescale_beta", "beta factor for rescaling diagonal", "", r_beta, &r_beta, 0);
  TaoOptionInt("-tao_lmm_scalar_history", "amount of history for scalar scaling", "", scalar_history, &scalar_history, 0);
  TaoOptionInt("-tao_lmm_rescale_history", "amount of history for rescaling diagonal", "", rescale_history, &rescale_history, 0);
  TaoOptionDouble("-tao_lmm_eps", "rejection tolerance", "", eps, &eps, 0);
  TaoOptionList("-tao_lmm_scale_type", "scale type", "", Scale_Table, Scale_Types, Scale_Table[scaleType], &scaleType, 0);
  TaoOptionList("-tao_lmm_rescale_type", "rescale type", "", Rescale_Table, Rescale_Types, Rescale_Table[rScaleType], &rScaleType, 0);
  TaoOptionList("-tao_lmm_limit_type", "limit type", "", Limit_Table, Limit_Types, Limit_Table[limitType], &limitType, 0);
  TaoOptionDouble("-tao_lmm_delta_min", "minimum delta value", "", delta_min, &delta_min, 0);
  TaoOptionDouble("-tao_lmm_delta_max", "maximum delta value", "", delta_max, &delta_max, 0);
  TaoOptionsTail();

  // Complete configuration
  rescale_history = TaoMin(rescale_history, lm);

  // Perform allocations
  tv->CloneVecs(lm+1, &S);
  tv->CloneVecs(lm+1, &Y);
  tv->Clone(&D);
  tv->Clone(&U);
  tv->Clone(&V);
  tv->Clone(&W);
  tv->Clone(&P);
  tv->Clone(&Q);

  rho = new double[lm+1];
  beta = new double[lm+1];

  yy_history = new double[TaoMax(scalar_history, 1)];
  ys_history = new double[TaoMax(scalar_history, 1)];
  ss_history = new double[TaoMax(scalar_history, 1)];

  yy_rhistory = new double[TaoMax(rescale_history, 1)];
  ys_rhistory = new double[TaoMax(rescale_history, 1)];
  ss_rhistory = new double[TaoMax(rescale_history, 1)];

  // Finish initializations
  lmnow = 0;
  iter = 0;
  updates = 0;
  rejects = 0;
  delta = 1.0;

  Gprev = 0;
  Xprev = 0;

  scale = 0;

  H0 = 0;
  H0default = TAO_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::TaoLMVMMat"
TaoLMVMMat::~TaoLMVMMat()
{
  int info;

  info = S[0]->DestroyVecs(lm+1, S);
  info = Y[0]->DestroyVecs(lm+1, Y);
  info = TaoVecDestroy(D);
  info = TaoVecDestroy(U);
  info = TaoVecDestroy(V);
  info = TaoVecDestroy(W);
  info = TaoVecDestroy(P);
  info = TaoVecDestroy(Q);

  delete[] rho;
  delete[] beta;
  if (info) rho=0;

  delete[] yy_history;
  delete[] ys_history;
  delete[] ss_history;

  delete[] yy_rhistory;
  delete[] ys_rhistory;
  delete[] ss_rhistory;
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::Refine"
int TaoLMVMMat::Refine(TaoLMVMMat *coarse, TaoMat *op, TaoVec *fineX, TaoVec *fineG)
{
  double rhotemp, rhotol;
  double y0temp, s0temp;
  int  info;
  TaoInt i;

  TaoFunctionBegin;

  info = Reset(); CHKERRQ(info);

  for (i = 0; i < coarse->lmnow; ++i) {
    // Refine S[i] and Y[i]
    info = op->Multiply(coarse->S[i], S[lmnow]); CHKERRQ(info);
    info = op->Multiply(coarse->Y[i], Y[lmnow]); CHKERRQ(info);

    // Check to see if refined vectors are fine
    info = Y[lmnow]->Dot(S[lmnow], &rhotemp); CHKERRQ(info);
    info = Y[lmnow]->Dot(Y[lmnow], &y0temp); CHKERRQ(info);

    rhotol = eps * y0temp;
    if (rhotemp > rhotol) {
      rho[lmnow] = 1.0 / rhotemp;

      // Add information to the scaling
      switch(scaleType) {
      case Scale_None:
        break;
	
      case Scale_Scalar:
        // Compute s^T s 
        info = S[lmnow]->Dot(S[lmnow], &s0temp); CHKERRQ(info);

        // Save information for scalar scaling
        yy_history[updates % scalar_history] = y0temp;
        ys_history[updates % scalar_history] = rhotemp;
        ss_history[updates % scalar_history] = s0temp;

	sigma = coarse->sigma;
	break;

      case Scale_Broyden:
	switch(rScaleType) {
	case Rescale_None:
	  break;
	  
        case Rescale_Scalar:
        case Rescale_GL:
          // Compute s^T s 
          info = S[lmnow]->Dot(S[lmnow], &s0temp); CHKERRQ(info);

          // Save information for special cases of scalar rescaling
          yy_rhistory[updates % rescale_history] = y0temp;
          ys_rhistory[updates % rescale_history] = rhotemp;
          ss_rhistory[updates % rescale_history] = s0temp;
	  break;
	}

	// Obtain diagonal scaling and ensure positive definite
	info = op->Multiply(coarse->D, D); CHKERRQ(info);
        info = D->AbsoluteValue(); CHKERRQ(info);
	break;
      }

      // Successful update
      ++updates;
      ++lmnow;
    }
    else {
      ++rejects;
    }
  }

  // Save the number of iterations and previous values for x and g
  // Need to use the actual value for the gradient here and not
  // the refined version of the coarse gradient in order for the
  // Wolfe conditions to guarantee the update is acceptable.
  iter = coarse->iter;
  info = Xprev->CopyFrom(fineX); CHKERRQ(info);
  info = Gprev->CopyFrom(fineG); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::Reset"
int TaoLMVMMat::Reset()
{
  TaoInt i;

  TaoFunctionBegin;
  Gprev = Y[lm];
  Xprev = S[lm];
  for (i = 0; i < lm; ++i) {
    rho[i] = 0.0;
  }
  rho[0] = 1.0;

  // Set the scaling and diagonal scaling matrix
  switch(scaleType) {
  case Scale_None:
    sigma = 1.0;
    break;

  case Scale_Scalar:
    sigma = delta;
    break;

  case Scale_Broyden:
    D->SetToConstant(delta);
    break;
  }

  iter = 0;
  updates = 0;
  lmnow = 0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::Presolve"
int TaoLMVMMat::Presolve()
{
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::Solve"
int TaoLMVMMat::Solve(TaoVec *tb, TaoVec *dx, TaoTruth *tt)
{
  double   sq, yq, dd;
  int      info;
  TaoInt ll;
  TaoTruth scaled;

  TaoFunctionBegin;
  if (lmnow < 1) {
    rho[0] = 1.0;
  }

  info = dx->CopyFrom(tb); CHKERRQ(info);
  for (ll = 0; ll < lmnow; ++ll) {
    info = dx->Dot(S[ll], &sq); CHKERRQ(info);
    beta[ll] = sq * rho[ll];
    info = dx->Axpy(-beta[ll], Y[ll]); CHKERRQ(info);
  }

  scaled = TAO_FALSE;
  if (!scaled && H0) {
    // use gradient as initial guess
    TaoVec *ApproxSoln; 
    info = tb->Clone(&ApproxSoln); CHKERRQ(info);
    info = ApproxSoln->CopyFrom(tb); CHKERRQ(info);
    info = H0->Solve(dx, ApproxSoln, tt); CHKERRQ(info);
    if (*tt) {
      info = dx->Dot(ApproxSoln, &dd); CHKERRQ(info);
      if ((dd > 0.0) && !TaoInfOrNaN(dd)) {
	// Accept Hessian solve
        PetscPrintf(PETSC_COMM_WORLD,"TaoLMVMMat: Accept Hessian Solve...\n");
	info = dx->CopyFrom(ApproxSoln); CHKERRQ(info);
	scaled = TAO_TRUE;
      }
      else{
        PetscPrintf(PETSC_COMM_WORLD,"TaoLMVMMat: Reject Hessian Solve, dd %22.12e \n\n", dd );
      }
    }
    info = TaoVecDestroy(ApproxSoln); CHKERRQ(info);
  }

  if (!scaled && scale) {
    info = U->PointwiseMultiply(dx, scale); CHKERRQ(info);
    info = dx->Dot(U, &dd); CHKERRQ(info);
    if ((dd > 0.0) && !TaoInfOrNaN(dd)) {
      // Accept scaling
      info = dx->CopyFrom(U); CHKERRQ(info);
      scaled = TAO_TRUE;
    }
  }
  
  if (!scaled) {
    switch(scaleType) {
    case Scale_None:
      break;

    case Scale_Scalar:
      info = dx->Scale(sigma); CHKERRQ(info);
      break;
  
    case Scale_Broyden:
      info = dx->PointwiseMultiply(dx, D); CHKERRQ(info);
      break;
    }
  } 

  for (ll = lmnow-1; ll >= 0; --ll) {
    info = dx->Dot(Y[ll], &yq); CHKERRQ(info);
    info = dx->Axpy(beta[ll] - yq * rho[ll], S[ll]); CHKERRQ(info);
  }

  *tt = TAO_TRUE;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::Update"
int TaoLMVMMat::Update(TaoVec *x, TaoVec *g)
{
  double rhotemp, rhotol;
  double y0temp, s0temp;
  double yDy, yDs, sDs;
  double sigmanew, denom;
  int  info;
  TaoInt i;

  double yy_sum, ys_sum, ss_sum;

  TaoFunctionBegin;

  if (0 == iter) {
    info = Reset(); CHKERRQ(info);
  } 
  else {
    info = Gprev->Aypx(-1.0, g); CHKERRQ(info);
    info = Xprev->Aypx(-1.0, x); CHKERRQ(info);

    info = Gprev->Dot(Xprev, &rhotemp); CHKERRQ(info);
    info = Gprev->Dot(Gprev, &y0temp); CHKERRQ(info);

    rhotol = eps * y0temp;
    if (rhotemp > rhotol) {
      ++updates;

      lmnow = TaoMin(lmnow+1, lm);
      for (i = lm-1; i >= 0; --i) {
	S[i+1] = S[i];
	Y[i+1] = Y[i];
	rho[i+1] = rho[i];
      }
      S[0] = Xprev;
      Y[0] = Gprev;
      rho[0] = 1.0 / rhotemp;

      // Compute the scaling
      switch(scaleType) {
      case Scale_None:
        break;

      case Scale_Scalar:
        // Compute s^T s 
        info = Xprev->Dot(Xprev, &s0temp); CHKERRQ(info);

	// Scalar is positive; safeguards are not required.

        // Save information for scalar scaling
        yy_history[(updates - 1) % scalar_history] = y0temp;
        ys_history[(updates - 1) % scalar_history] = rhotemp;
        ss_history[(updates - 1) % scalar_history] = s0temp;

        // Compute summations for scalar scaling
        yy_sum = 0;	// No safeguard required; y^T y > 0
        ys_sum = 0;	// No safeguard required; y^T s > 0
        ss_sum = 0;	// No safeguard required; s^T s > 0
        for (i = 0; i < TaoMin(updates, scalar_history); ++i) {
          yy_sum += yy_history[i];
          ys_sum += ys_history[i];
          ss_sum += ss_history[i];
        }

        if (0.0 == s_alpha) {
	  // Safeguard ys_sum 
	  if (0.0 == ys_sum) {
            ys_sum = TAO_ZER_SAFEGUARD;
          }

          sigmanew = ss_sum / ys_sum;
        }
        else if (0.5 == s_alpha) {
	  // Safeguard yy_sum 
	  if (0.0 == yy_sum) {
            yy_sum = TAO_ZER_SAFEGUARD;
          }

          sigmanew = sqrt(ss_sum / yy_sum);
        }
        else if (1.0 == s_alpha) {
	  // Safeguard yy_sum 
	  if (0.0 == yy_sum) {
            yy_sum = TAO_ZER_SAFEGUARD;
          }

          sigmanew = ys_sum / yy_sum;
        }
        else {
	  denom = 2*s_alpha*yy_sum;

          // Safeguard denom
	  if (0.0 == denom) {
            denom = TAO_ZER_SAFEGUARD;
          }

          sigmanew = ((2*s_alpha-1)*ys_sum + 
                      sqrt((2*s_alpha-1)*(2*s_alpha-1)*ys_sum*ys_sum - 
                           4*s_alpha*(s_alpha-1)*yy_sum*ss_sum)) / denom;
        }

	switch(limitType) {
	case Limit_Average:
          if (1.0 == mu) {
            sigma = sigmanew;
          }
          else if (mu) {
            sigma = mu * sigmanew + (1.0 - mu) * sigma;
          }
	  break;

        case Limit_Relative:
          if (mu) {
            sigma = TaoMid((1.0 - mu) * sigma, sigmanew, (1.0 + mu) * sigma);
          }
          break;

        case Limit_Absolute:
          if (nu) {
            sigma = TaoMid(sigma - nu, sigmanew, sigma + nu);
          }
          break;

        default:
	  sigma = sigmanew;
	  break;
        }
        break;

      case Scale_Broyden:
        // Original version
        // Combine DFP and BFGS

	// This code appears to be numerically unstable.  We use the
	// original version because this was used to generate all of
	// the data and because it may be the least unstable of the
	// bunch.

        // P = Q = inv(D);
        info = P->CopyFrom(D); CHKERRQ(info);
        info = P->Reciprocal(); CHKERRQ(info);
        info = Q->CopyFrom(P); CHKERRQ(info);

        // V = y*y
        info = V->PointwiseMultiply(Gprev, Gprev); CHKERRQ(info);

        // W = inv(D)*s
        info = W->PointwiseMultiply(Xprev, P); CHKERRQ(info);
        info = W->Dot(Xprev, &sDs); CHKERRQ(info);

        // Safeguard rhotemp and sDs
        if (0.0 == rhotemp) {
          rhotemp = TAO_ZER_SAFEGUARD;
        }

        if (0.0 == sDs) {
          sDs = TAO_ZER_SAFEGUARD;
        }

        if (1.0 != phi) {
          // BFGS portion of the update
          // U = (inv(D)*s)*(inv(D)*s)
          info = U->PointwiseMultiply(W, W); CHKERRQ(info);

          // Assemble
          info = P->Axpy(1.0 / rhotemp, V); CHKERRQ(info);
          info = P->Axpy(-1.0 / sDs, U); CHKERRQ(info);
        }

        if (0.0 != phi) {
          // DFP portion of the update
          // U = inv(D)*s*y
          info = U->PointwiseMultiply(W, Gprev); CHKERRQ(info);

          // Assemble
          info = Q->Axpy(1.0 / rhotemp + sDs / (rhotemp * rhotemp), V); CHKERRQ(info);
          info = Q->Axpy(-2.0 / rhotemp, U); CHKERRQ(info);
        }

        if (0.0 == phi) {
          info = U->CopyFrom(P); CHKERRQ(info);
        }
        else if (1.0 == phi) {
          info = U->CopyFrom(Q); CHKERRQ(info);
        }
        else {
          // Broyden update
          info = U->Waxpby(1.0 - phi, P, phi, Q); CHKERRQ(info);
        }

	// Obtain inverse and ensure positive definite
	info = U->Reciprocal(); CHKERRQ(info);
        info = U->AbsoluteValue(); CHKERRQ(info);

        // Checking the diagonal scaling for not a number and infinity
	// should not be necessary for the Broyden update
        // info = U->Dot(U, &sDs); CHKERRQ(info);
        // if (sDs != sDs) {
        //   // not a number
        //   info = U->SetToConstant(TAO_ZER_SAFEGUARD); CHKERRQ(info);
        // }
        // else if ((sDs - sDs) != 0.0)) {
        //   // infinity
        //   info = U->SetToConstant(TAO_INF_SAFEGUARD); CHKERRQ(info);
        // }

	switch(rScaleType) {
	case Rescale_None:
	  break;

        case Rescale_Scalar:
        case Rescale_GL:
	  if (rScaleType == Rescale_GL) {
	    // Gilbert and Lemarachal use the old diagonal
	    info = P->CopyFrom(D); CHKERRQ(info);
          }
          else {
	    // The default version uses the current diagonal
	    info = P->CopyFrom(U); CHKERRQ(info);
          }

          // Compute s^T s 
          info = Xprev->Dot(Xprev, &s0temp); CHKERRQ(info);

          // Save information for special cases of scalar rescaling
          yy_rhistory[(updates - 1) % rescale_history] = y0temp;
          ys_rhistory[(updates - 1) % rescale_history] = rhotemp;
          ss_rhistory[(updates - 1) % rescale_history] = s0temp;

          if (0.0 == r_beta) {
            if (1 == TaoMin(updates, rescale_history)) {
              // Compute summations for scalar scaling
	      info = W->PointwiseDivide(S[0], P); CHKERRQ(info);
  
              info = W->Dot(Y[0], &ys_sum); CHKERRQ(info);
              info = W->Dot(W, &ss_sum); CHKERRQ(info);
  
              yy_sum += yy_rhistory[0];
            }
            else {
              info = Q->CopyFrom(P); CHKERRQ(info);
              info = Q->Reciprocal(); CHKERRQ(info);

              // Compute summations for scalar scaling
              yy_sum = 0;	// No safeguard required
              ys_sum = 0;	// No safeguard required
              ss_sum = 0;	// No safeguard required
              for (i = 0; i < TaoMin(updates, rescale_history); ++i) {
	        info = W->PointwiseMultiply(S[i], Q); CHKERRQ(info);
                info = W->Dot(Y[i], &yDs); CHKERRQ(info);
                ys_sum += yDs;

                info = W->Dot(W, &sDs); CHKERRQ(info);
                ss_sum += sDs;
  
                yy_sum += yy_rhistory[i];
              }
            }
          }
          else if (0.5 == r_beta) {
            if (1 == TaoMin(updates, rescale_history)) {
	      info = V->PointwiseMultiply(Y[0], P); CHKERRQ(info);
	      info = V->Dot(Y[0], &yy_sum); CHKERRQ(info);

	      info = W->PointwiseDivide(S[0], P); CHKERRQ(info);
              info = W->Dot(S[0], &ss_sum); CHKERRQ(info);

              ys_sum = ys_rhistory[0];
            }
            else {
              info = Q->CopyFrom(P); CHKERRQ(info);
              info = Q->Reciprocal(); CHKERRQ(info);

              // Compute summations for scalar scaling
              yy_sum = 0;	// No safeguard required
              ys_sum = 0;	// No safeguard required
              ss_sum = 0;	// No safeguard required
              for (i = 0; i < TaoMin(updates, rescale_history); ++i) {
	        info = V->PointwiseMultiply(Y[i], P); CHKERRQ(info);
	        info = V->Dot(Y[i], &yDy); CHKERRQ(info);
                yy_sum += yDy;

	        info = W->PointwiseMultiply(S[i], Q); CHKERRQ(info);
                info = W->Dot(S[i], &sDs); CHKERRQ(info);
                ss_sum += sDs;

                ys_sum += ys_rhistory[i];
              }
            }
	  }
          else if (1.0 == r_beta) {
            // Compute summations for scalar scaling
            yy_sum = 0;	// No safeguard required
            ys_sum = 0;	// No safeguard required
            ss_sum = 0;	// No safeguard required
            for (i = 0; i < TaoMin(updates, rescale_history); ++i) {
	      info = V->PointwiseMultiply(Y[i], P); CHKERRQ(info);
              info = V->Dot(S[i], &yDs); CHKERRQ(info);
              ys_sum += yDs;

	      info = V->Dot(V, &yDy); CHKERRQ(info);
              yy_sum += yDy;

              ss_sum += ss_rhistory[i];
            }
          }
          else {
	    info = Q->CopyFrom(P); CHKERRQ(info);

            info = P->Pow(r_beta); CHKERRQ(info);
	    info = Q->PointwiseDivide(P, Q); CHKERRQ(info);

            // Compute summations for scalar scaling
            yy_sum = 0;	// No safeguard required
            ys_sum = 0;	// No safeguard required
            ss_sum = 0;	// No safeguard required
            for (i = 0; i < TaoMin(updates, rescale_history); ++i) {
	      info = V->PointwiseMultiply(P, Y[i]); CHKERRQ(info);
	      info = W->PointwiseMultiply(Q, S[i]); CHKERRQ(info);

	      info = V->Dot(V, &yDy); CHKERRQ(info);
              info = V->Dot(W, &yDs); CHKERRQ(info);
              info = W->Dot(W, &sDs); CHKERRQ(info);

              yy_sum += yDy;
              ys_sum += yDs;
              ss_sum += sDs;
            }
          }

          if (0.0 == r_alpha) {
	    // Safeguard ys_sum 
	    if (0.0 == ys_sum) {
              ys_sum = TAO_ZER_SAFEGUARD;
            }

            sigmanew = ss_sum / ys_sum;
          }
          else if (0.5 == r_alpha) {
	    // Safeguard yy_sum 
	    if (0.0 == yy_sum) {
              yy_sum = TAO_ZER_SAFEGUARD;
            }

            sigmanew = sqrt(ss_sum / yy_sum);
          }
          else if (1.0 == r_alpha) {
	    // Safeguard yy_sum 
	    if (0.0 == yy_sum) {
              yy_sum = TAO_ZER_SAFEGUARD;
            }

            sigmanew = ys_sum / yy_sum;
          }
          else {
	    denom = 2*r_alpha*yy_sum;

            // Safeguard denom
	    if (0.0 == denom) {
              denom = TAO_ZER_SAFEGUARD;
            }

            sigmanew = ((2*r_alpha-1)*ys_sum +
                        sqrt((2*r_alpha-1)*(2*r_alpha-1)*ys_sum*ys_sum -
                             4*r_alpha*(r_alpha-1)*yy_sum*ss_sum)) / denom;
          }

	  // If Q has small values, then Q^(r_beta - 1)
          // can have very large values.  Hence, ys_sum
          // and ss_sum can be infinity.  In this case,
	  // sigmanew can either be not-a-number or infinity.

          if (TaoInfOrNaN(sigmanew)) {
            // sigmanew is not-a-number; skip rescaling
          }
          else if (!sigmanew) {
	    // sigmanew is zero; this is a bad case; skip rescaling
          }
          else {
	    // sigmanew is positive
            info = U->Scale(sigmanew); CHKERRQ(info);
          }
	  break;
	}

        // Modify for previous information
	switch(limitType) {
	case Limit_Average:
	  if (1.0 == mu) {
            info = D->CopyFrom(U); CHKERRQ(info);
          }
          else if (mu) {
            info = D->Axpby(mu, U, 1.0 - mu); CHKERRQ(info);
          }
	  break;
 
        case Limit_Relative:
	  if (mu) {
	    info = P->ScaleCopyFrom(1.0 - mu, D); CHKERRQ(info);
	    info = Q->ScaleCopyFrom(1.0 + mu, D); CHKERRQ(info);
            info = D->Median(P, U, Q); CHKERRQ(info);
          }
          break;

        case Limit_Absolute:
	  if (nu) {
	    info = P->CopyFrom(D); CHKERRQ(info);
	    info = P->AddConstant(-nu); CHKERRQ(info);
   	    info = Q->CopyFrom(D); CHKERRQ(info);
	    info = Q->AddConstant(nu); CHKERRQ(info);
            info = D->Median(P, U, Q); CHKERRQ(info);
          }
	  break;

        default:
          info = D->CopyFrom(U); CHKERRQ(info);
	  break;
        } 
	break;
      }

      Xprev = S[lm]; Gprev = Y[lm];
    } 
    else { 
      ++rejects;
    }
  }
  
  ++iter;
  info = Xprev->CopyFrom(x); CHKERRQ(info);
  info = Gprev->CopyFrom(g); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::View"
int TaoLMVMMat::View()
{
  TaoFunctionBegin;
  PetscErrorCode info;
  info=PetscPrintf(PETSC_COMM_WORLD,"TaoLMVMMat: lm %d, lmnow %d, iter %d, updates %d, rejects %d\n",lm, lmnow, iter, updates, rejects);
  std::cout << std::flush ;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::SetDelta()"
int TaoLMVMMat::SetDelta(double d)
{
  TaoFunctionBegin;
  delta = TaoAbsDouble(d);
  delta = TaoMax(delta_min, delta);
  delta = TaoMin(delta_max, delta);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::SetDelta()"
int TaoLMVMMat::SetScale(TaoVec *s)
{
  TaoFunctionBegin;
  scale = s;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::SetH0()"
int TaoLMVMMat::SetH0(TaoMat *HH0)
{
  TaoFunctionBegin;
  if (H0) { 
    H0 = 0;
  }
  H0default = TAO_TRUE;

  if (HH0) {
    H0 = HH0;
    H0default = TAO_FALSE;
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::GetX0"
int TaoLMVMMat::GetX0(TaoVec *x)
{
  int i,info;

  TaoFunctionBegin;
  info = x->CopyFrom(Xprev); CHKERRQ(info);
  for (i = 0; i < lmnow; ++i) {
    info = x->Axpy(-1.0, S[i]); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLMVMMat::InitialApproximation"
int TaoLMVMMat::InitialApproximation(TaoVec *x)
{
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

