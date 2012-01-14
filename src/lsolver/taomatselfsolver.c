#include "tao_general.h"    /*I "tao_solver.h"  I*/
#include "taomatselfsolver.h"
#include "taomat.h"

#undef __FUNCT__
#define __FUNCT__ "TaoMatSelfSolver::PreSolve"
int TaoMatSelfSolver::PreSolve(TaoMat* M){
  int info;
  TaoFunctionBegin;
  this->tmoperator=M;
  info=M->Presolve(); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatSelfSolver::Solve"
int TaoMatSelfSolver::Solve(TaoVec* b, TaoVec* x, TaoTruth *flag){
  int info;
  TaoFunctionBegin;
  if (!this->tmoperator){
    SETERRQ(56,"No PreSolve() operation called or invalide matrix.");
  }
  info = this->tmoperator->Solve(b,x,flag);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatSelfSolver::GetNumberIterations"
int TaoMatSelfSolver::GetNumberIterations(int * iters){
  TaoFunctionBegin;
  *iters=1;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatSelfSolver::SetTolerances"
int TaoMatSelfSolver::SetTolerances(double rtol, double atol, double dtol, int maxits){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatSelfSolver::SetOptions"
int TaoMatSelfSolver::SetOptions(){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatSelfSolver::View"
int TaoMatSelfSolver::View(){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

