/*$Id$*/

/*
 Context for limited memory variable metric method for bound constrained 
 optimization.
*/

#ifndef __TAO_BLMVM_H
#define __TAO_BLMVM_H
#include "tao_solver.h"
#include "src/matrix/lmvmmat.h"

typedef struct {

  TaoLMVMMat *M;

  TaoVec *G;
  TaoVec *GP;
  TaoVec *D;
  TaoVec *XL;
  TaoVec *XU;

  TaoInt n_free;
  TaoInt n_bind;

  TaoInt grad;
  TaoInt reset;
} TAO_BLMVM;

#endif
