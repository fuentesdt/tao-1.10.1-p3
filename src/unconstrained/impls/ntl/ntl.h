/*$Id: ntl.h,v 1.2 2008-08-18 19:50:51 sarich Exp $*/

// Context for a Newton trust-region, line-search method for unconstrained 
// minimization

#ifndef __TAO_NTL_H
#define __TAO_NTL_H
#include "tao_solver.h"
#include "src/matrix/lmvmmat.h"

typedef struct {
  TaoLMVMMat *M;

  TaoVec *G;
  TaoVec *D;
  TaoVec *W;

  TaoVec *Xold;
  TaoVec *Gold;
  TaoVec *Diag;

  // Parameters when updating the trust-region radius based on steplength
  double nu1;		// used to compute trust-region radius
  double nu2;		// used to compute trust-region radius
  double nu3;		// used to compute trust-region radius
  double nu4;		// used to compute trust-region radius

  double omega1;        // factor used for trust-region update
  double omega2;        // factor used for trust-region update
  double omega3;        // factor used for trust-region update
  double omega4;        // factor used for trust-region update
  double omega5;        // factor used for trust-region update

  // if   step < nu1  		(very bad step)
  //   radius = omega1 * min(norm(d), radius)
  // elif step < nu2		(bad step)
  //   radius = omega2 * min(norm(d), radius)
  // elif step < nu3		(okay step)
  //   radius = omega3 * radius;
  // elif step < nu4		(good step)
  //   radius = max(omega4 * norm(d), radius)
  // else 			(very good step)
  //   radius = max(omega5 * norm(d), radius)
  // fi
 
  // Parameters when updating the trust-region radius based on reduction
  double eta1;		// used to compute trust-region radius
  double eta2;		// used to compute trust-region radius
  double eta3;		// used to compute trust-region radius
  double eta4;		// used to compute trust-region radius

  double alpha1;        // factor used for trust-region update
  double alpha2;        // factor used for trust-region update
  double alpha3;        // factor used for trust-region update
  double alpha4;        // factor used for trust-region update
  double alpha5;        // factor used for trust-region update

  // kappa = ared / pred
  // if   kappa < eta1 		(very bad step)
  //   radius = alpha1 * min(norm(d), radius)
  // elif kappa < eta2		(bad step)
  //   radius = alpha2 * min(norm(d), radius)
  // elif kappa < eta3		(okay step)
  //   radius = alpha3 * radius;
  // elif kappa < eta4		(good step)
  //   radius = max(alpha4 * norm(d), radius)
  // else 			(very good step)
  //   radius = max(alpha5 * norm(d), radius)
  // fi
 
  // Parameters when updating the trust-region radius based on interpolation
  double mu1;		// used for model agreement in interpolation
  double mu2;		// used for model agreement in interpolation

  double gamma1;	// factor used for interpolation
  double gamma2;	// factor used for interpolation
  double gamma3;	// factor used for interpolation
  double gamma4;	// factor used for interpolation

  double theta;		// factor used for interpolation

  // kappa = ared / pred
  // if   kappa >= 1.0 - mu1	(very good step)
  //   choose tau in [gamma3, gamma4]
  //   radius = max(tau * norm(d), radius)
  // elif kappa >= 1.0 - mu2    (good step)
  //   choose tau in [gamma2, gamma3]
  //   if (tau >= 1.0)
  //     radius = max(tau * norm(d), radius)
  //   else
  //     radius = tau * min(norm(d), radius)
  //   fi
  // else 			(bad step)
  //   choose tau in [gamma1, 1.0]
  //   radius = tau * min(norm(d), radius)
  // fi
 
  // Parameters when initializing trust-region radius based on interpolation
  double mu1_i;		// used for model agreement in interpolation
  double mu2_i;		// used for model agreement in interpolation

  double gamma1_i;	// factor used for interpolation
  double gamma2_i;	// factor used for interpolation
  double gamma3_i;	// factor used for interpolation
  double gamma4_i;	// factor used for interpolation

  double theta_i;	// factor used for interpolation

  // Other parameters
  double min_radius;    // lower bound on initial radius value
  double max_radius;    // upper bound on trust region radius
  double epsilon;       // tolerance used when computing ared/pred

  TaoInt trust;		// Trust-region steps accepted
  TaoInt newt;		// Newton directions attempted
  TaoInt bfgs;		// BFGS directions attempted
  TaoInt sgrad;		// Scaled gradient directions attempted
  TaoInt grad;		// Gradient directions attempted

  TaoInt ksp_type;		// KSP method for the code
  TaoInt pc_type;		// Preconditioner for the code
  TaoInt bfgs_scale_type;	// Scaling matrix to used for the bfgs preconditioner
  TaoInt init_type;	// Trust-region initialization method
  TaoInt update_type;      // Trust-region update method
} TAO_NTL;

#endif
