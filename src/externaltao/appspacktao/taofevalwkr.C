/* taofevalwkr.C */

#include "fevalwkr.H"
#include "taofevalwkr.h"
#include "tao_solver.h"



// *** This MUST be called 'myfevalwkr' ***
//==========================================
#undef __FUNCT__
#define __FUNCT__ "myfevalwkr"
FevalWkr* myfevalwkr()
{
  return new TaoFevalWkr;
}


//==========================================
#undef __FUNCT__
#define __FUNCT__ "TaoFevalWkr::doFeval"
//double TaoFevalWkr::doFeval(double * truex, bool &isf)
double TaoFevalWkr::doFeval(vector<double>& truex, bool &isf)
{
  double fval;
  int info;
  int ndim,i;
  double *x;
  TaoVec *XX;
  static int iters = 0;
  int n;
  TaoTerminateReason reason;

  // tao, appsPtr are class members

  // Convert the array truex to a TaoVec 
  info = TaoGetSolution(tao,&XX); CHKERRQ(info);
  info = XX->GetArray(&x,&n); CHKERRQ(info);
  for (i=0; i<n;i++)
    x[i] = truex[i];
  info = XX->RestoreArray(&x,&n); CHKERRQ(info);


  /* Now get the function value */
  info = TaoComputeFunction(tao, XX, &fval); CHKERRQ(info);
  appsPtr->fval = fval;

  iters++;


  // update tao;

  info = TaoMonitor(tao, iters, fval, 1, 0, 1, &reason); CHKERRQ(info);
  isf = true; // Not sure what this is for
  
  
  return fval;
}
  
  
