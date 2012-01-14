#include "tao_general.h"
#include "taomat_ga.h"
#include "taovec_ga.h"
#include "taolinearsolver_ga.h"

TaoLinearSolverGa::TaoLinearSolverGa():TaoLinearSolver() {
  linear_iters=0;
  mm = 0;
}

int 
TaoLinearSolverGa::SetOperator(TaoMat* m1){
  TaoFunctionBegin;
  mm = ((TaoMatGa *)m1)->GetMat();
  TaoFunctionReturn(0);
}


int 
TaoLinearSolverGa::SetTrustRadius(double rad){
  TaoFunctionBegin;
  GA_Error ((char*)"TaoLinearSolverGA::SetTrustRadius() not defined: ", 0L);
  TaoFunctionReturn(1);
}

int 
TaoLinearSolverGa::Solve(TaoVec* v1, TaoVec* v2){
  GAVec g_vv=((TaoVecGa *)v1)->GetVec();
  GAVec g_ww=((TaoVecGa *)v2)->GetVec();
  TaoFunctionBegin;
  GA_Copy(g_vv, g_ww);
  GA_Lu_solve('N', mm, g_ww);
  TaoFunctionReturn(0);
}

int 
TaoLinearSolverGa::GetNumberIterations(int * iters){
  TaoFunctionBegin;
  *iters=linear_iters;
  TaoFunctionReturn(0);
}

int 
TaoLinearSolverGa::SetTolerances(double rtol, double atol, double dtol, 
				 int maxits){
  TaoFunctionBegin;
  GA_Error ((char*)"TaoLinearSolverGA::SetTolerances() not defined: ", 0L);
  TaoFunctionReturn(1);
}

int 
TaoLinearSolverGa::View(){
  TaoFunctionBegin;
  GA_Error ((char*)"TaoLinearSolverGA::View() not defined: ", 0L);
  TaoFunctionReturn(1);
}

int 
TaoLinearSolverGa::SetOptions(){
  TaoFunctionBegin;
  GA_Error ((char*)"TaoLinearSolverGA::SetOptions() not defined: ", 0L);
  TaoFunctionReturn(1);
}

int 
TaoLinearSolverGa::Duplicate(TaoLinearSolver**M){
  TaoFunctionBegin;
  *M = new TaoLinearSolverGa(*this);
  TaoFunctionReturn(0);
}
