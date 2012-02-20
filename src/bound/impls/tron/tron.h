/*$Id$*/

#ifndef __TAO_TRON_H
#define __TAO_TRON_H

#include "tao_solver.h"
#include "src/matrix/lmvmmat.h"

typedef struct {

  /* Parameters */
  double pg_ftol;
  double actred;
  double f_new;
 
  double eta1,eta2,eta3,eta4;
  double sigma1,sigma2,sigma3;

  TaoInt maxgpits;

  /* Problem variables, vectors and index sets */
  double stepsize;
  double pgstepsize;

  /* Problem statistics */

  TaoInt n;   /* Dimension of the Problem */
  double delta;  /* Trust region size */
  double gnorm;
  double f;

  TaoInt total_cgits;
  TaoInt cg_iterates;
  TaoInt total_gp_its;
  TaoInt gp_iterates;
  TaoInt cgits;

  TaoVec* DXFree;
  TaoVec* R;

  TaoVec* X;
  TaoVec* G;
  TaoVec* PG;

  TaoVec* DX;
  TaoVec* X_New;
  TaoVec* G_New;
  TaoVec* XU;
  TaoVec* XL;
  TaoVec* Work;
  
  TaoMat* Hsub;
  TaoMat* H;
  TaoLMVMMat *M;

  TaoIndexSet *TT;
  TaoIndexSet *Free_Local;  /* Indices of local variables equal to lower bound */
  TaoIndexSet *Lower_Local;  /* Indices of local variables equal to lower bound */
  TaoIndexSet *Upper_Local;  /* Indices of local variables equal to lower bound */

  TaoInt n_free;       /* Number of free variables */
  TaoInt n_upper;
  TaoInt n_lower;
  TaoInt n_bind;       /* Number of binding varibles */

} TAO_TRON;

#endif

