#ifndef TAOABCFLOATAPP_H
#define TAOABCFLOATAPP_H

#include "taoappobject.h"
#include "src/vector/tvecsingle.h"

class TaoABCFloatApplication: public TaoApplication{

 private:

 public:
  TaoABCFloatApplication();
  ~TaoABCFloatApplication();

  TaoVecFloatArray *taox;   /* The Variable Vector */
  int SetNumberOfVariables(TaoInt);
  int GetSolutionAndGradient(float*&, float*&, TaoInt&);

  /* Function and Gradient */
  int GetVariableVector(TaoVec **);

  int EvaluateObjectiveFunction(TaoVec *, double *);
  int EvaluateGradient(TaoVec *, TaoVec *);
  int EvaluateObjectiveAndGradient(TaoVec *xx, double *ff, TaoVec *gg);

  virtual int  ComputeObjectiveAndGradient(float*,TaoInt,float*,float *)=0;

  /* Set Variable Bounds */
  int EvaluateVariableBounds(TaoVec *xxll, TaoVec *xxuu);
  virtual int ComputeVariableBounds(float *xl, float *xu, TaoInt n);

  /* Initialize Starting Point */
  int InitializeVariables(TaoVec *x);
  virtual int StartingPoint(float *x, TaoInt n);

};

#endif


