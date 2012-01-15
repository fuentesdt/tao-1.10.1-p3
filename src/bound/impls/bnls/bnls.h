/*$Id: s.tron.h 1.35 02/08/16 16:55:59-05:00 benson@rockies.mcs.anl.gov $*/

#ifndef __TAO_BNLS_H
#define __TAO_BNLS_H
#include "tao_solver.h"
#include "src/matrix/lmvmmat.h"

typedef struct{

  /* Problem statistics */

  double gamma;		     /* damping parameter */
  double gamma_factor;	     /* damping parameter */

  TaoVec* DXFree;
  TaoVec* R;

  TaoVec* X;
  TaoVec* G, *PG;

  TaoVec* DX;
  TaoVec* XU;
  TaoVec* XL;
  TaoVec* Work;
  
  TaoMat* Hsub;
  TaoMat* H;
  TaoLMVMMat *M;

  TaoIndexSet *FreeVariables;  /* Indices of local variables equal to lower bound */


}TAO_BNLS;

#endif






