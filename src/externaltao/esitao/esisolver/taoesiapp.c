#include "../interface/taovec_esi.h"
#include "../interface/taomat_esi.h"
#include "../interface/taolinearsolver_esi.h"
#include "esiapplication.h"
#include "taoesiapp.h"
#include "tao_solver.h"
#include "src/tao_impl.h"
#include <stdlib.h>


#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::TaoESIApplication"
TaoESIApplication::TaoESIApplication(esi::ESIApplication* eapp){
  esiapp=eapp;
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::~TaoESIApplication"
TaoESIApplication::~TaoESIApplication(){
  return;
}


#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::EvaluateObjectiveFunction"
int TaoESIApplication::EvaluateObjectiveFunction(TaoVec *xx, double *ff){
  int info;
  esi::Vector<double,int> *x=(esi::Vector<double,int>*)(xx->VecObject);
  TaoFunctionBegin;
  info=esiapp->evaluateObjectiveFunction(x,ff); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::EvaluateGradient"
int TaoESIApplication::EvaluateGradient(TaoVec *xx, TaoVec *gg){
  int info;
  esi::Vector<double,int> *x=(esi::Vector<double,int>*)(xx->VecObject);
  esi::Vector<double,int> *g=(esi::Vector<double,int>*)(gg->VecObject);
  TaoFunctionBegin;
  info=esiapp->evaluateGradient(x,g); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::EvaluateObjectiveAndGradient"
int TaoESIApplication::EvaluateObjectiveAndGradient(TaoVec *xx, double *ff, TaoVec *gg){
  int info;
  esi::Vector<double,int> *x=(esi::Vector<double,int>*)(xx->VecObject);
  esi::Vector<double,int> *g=(esi::Vector<double,int>*)(gg->VecObject);
  TaoFunctionBegin;
  info=esiapp->evaluateObjectiveAndGradient(x,ff,g);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::InitializeVariables"
int TaoESIApplication::InitializeVariables(TaoVec *x){
  TaoFunctionBegin;
  //  info=x->setToZero();CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::GetVariableVector"
int TaoESIApplication::GetVariableVector(TaoVec **xx){
  int info;
  esi::Vector<double,int> *x;
  TaoFunctionBegin;
  info=esiapp->getVariableVector(&x); CHKERRQ(info);
  *xx=new TaoVecESI(x);
  //  info=x->deleteReference();
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::GetVariableBounds"
int TaoESIApplication::EvaluateVariableBounds(TaoVec *xxll, TaoVec *xxuu){
  int info;
  esi::Vector<double,int> *xl=(esi::Vector<double,int>*)(xxll->VecObject);
  esi::Vector<double,int> *xu=(esi::Vector<double,int>*)(xxuu->VecObject);

  TaoFunctionBegin;
  info=esiapp->evaluateBounds(xl,xu); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::GetHessianMatrix"
int TaoESIApplication::GetHessianMatrix(TaoMat **HH){
  int info;
  esi::Operator<double,int> *H;
  TaoFunctionBegin;
  info=esiapp->getHessianMatrix(&H); CHKERRQ(info);
  *HH=0;
  *HH=new TaoMatESI(H);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::EvaluateHessian"
int TaoESIApplication::EvaluateHessian(TaoVec *xx, TaoMat *HH){
  int info;
  TaoFunctionBegin;
  esi::Vector<double,int> *x=(esi::Vector<double,int>*)(xx->VecObject);
  esi::Operator<double,int> *H=((TaoMatESI*)HH)->esiobj;
  info=esiapp->evaluateHessian(x,H); CHKERRQ(info);
  TaoFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "TaoESIApplication::Monitor"
int TaoESIApplication::Monitor(){
  int info;
  TaoFunctionBegin;
  info=esiapp->monitor(); CHKERRQ(info);
  TaoFunctionReturn(0);
}
