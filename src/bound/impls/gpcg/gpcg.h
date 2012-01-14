/*$Id$*/

#ifndef __TAO_GPCG_H
#define __TAO_GPCG_H
#include "tao_solver.h"
#include "src/bound/impls/gpcg/gpcglinesearch.h"

typedef struct{

  /* Parameters */
  double pg_ftol;
  double actred;
  double f_new;
  double minstep;
  double stepsize;
  double gnorm;

  double sigma1,sigma2,sigma3;

  TaoInt maxgpits;

  /* Problem variables, vectors and index sets */

  /* Problem statistics */

  TaoInt n;   /* Dimension of the Problem */

  TaoInt total_cgits;
  TaoInt cg_iterates;
  TaoInt total_gp_its;
  TaoInt gp_iterates;
  TaoInt cgits;

  TaoVec * G_New;
  TaoVec * DXFree;
  TaoVec * R;
  TaoVec * DX;
  TaoVec * X;
  TaoVec * X_New;
  TaoVec * G, *PG;
  TaoVec * XU;
  TaoVec * XL;
  TaoVec * Work;

  TaoMat *H;
  TaoVec * B;
  double c;
  
  double f;
  double step;
  TaoMat *Hsub;

  TaoIndexSet * Free_Local;  /* Indices of local variables equal to lower bound */
  TaoIndexSet * TT;  /* Indices of local variables equal to upper bound */

  TaoInt n_free;       /* Number of free variables */
  TaoInt n_upper;
  TaoInt n_lower;
  TaoInt n_bind;       /* Number of binding varibles */

}TAO_GPCG;

/* GPCG Routines */
int TaoGPCGComputeFunctionGradient(TAO_SOLVER, TaoVec *, double *, TaoVec *);

#endif






