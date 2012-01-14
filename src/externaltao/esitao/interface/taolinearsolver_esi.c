
#include "esi/ESI.h"
#include "taolinearsolver_esi.h"
#include "tao_general.h"

static int Tao_LinearSolve=0;
static int TAO_LINEARSOLVER_COOKIE=0;

TaoLinearSolverESI::TaoLinearSolverESI(esi::Solver<double,int>* S):TaoLinearSolver(){
  S->addReference();
  ksp=S;
  its=0;
  if (TAO_LINEARSOLVER_COOKIE==0){
    PetscLogClassRegister(&TAO_LINEARSOLVER_COOKIE,"TAO Application"); 
  }
  PetscLogEventRegister(&Tao_LinearSolve,"TaoLinearSolve",TAO_LINEARSOLVER_COOKIE);

  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoWrapESIKSP"
/*@C
   TaoWrapKSP - Create a new TaoLinearSolver object using PETSc KSP.

   Input Parameter:
+  S -  an ESI Linear solver
-  flg - set to TAO_TRUE if the KSP S should be destroyed when this 
   TaoLinearSolver is destroyed, and set to TAO_FALSE otherwise.

   Output Parameter:
.  SS - new TaoMat

   Level: beginner

@*/
int TaoWrapESIKSP( esi::Solver<double,int>* S, TaoLinearSolverESI ** SS){
  TaoFunctionBegin;
  *SS = new  TaoLinearSolverESI(S);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverGetKSP"
/*@C
   TaoLinearSolverGetKSP - If the TaoLinearSolver is of the TaoLinearSolverESI type, it gets the underlying PETSc KSP 

   Input Parameter:
.  TS - the TaoLinearSolver

   Output Parameter:
.  S -  the underlying PETSc KSP object

Note:
The function TaoMatESI::GetKSP() will also return the KSP S.

   Level: beginner

@*/
int TaoLinearSolverGetESIKSP( TaoLinearSolver *TS,  esi::Solver<double,int>** S){
  TaoFunctionBegin;
  *S=((TaoLinearSolverESI *)TS)->ksp;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverESI::PreSolve"
int TaoLinearSolverESI::PreSolve(TaoMat* m1)
{
  esi::Operator<double,int> *mm=(esi::Operator<double,int>*)(m1->VecObject);
  TaoFunctionBegin;
  info=mm->esiobj->setup(); CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverESI::Solve"
int TaoLinearSolverESI::Solve(TaoVec* tv, TaoVec* tw){
  esi::Vector<double,int> *vv=(esi::Vector<double,int>*)(tv->VecObject);
  esi::Vector<double,int> *ww=(esi::Vector<double,int>*)(tw->VecObject);
  esi::ErrorCode info;

  TaoFunctionBegin;
  its=1;
  info=ksp->solve(*vv,*ww); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverESI::GetNumberIterations"
int TaoLinearSolverESI::GetNumberIterations(int * iters){
  esi::SolverIterative<double,int> *SS;
  esi::ErrorCode info;
  TaoFunctionBegin;
  SS= dynamic_cast<esi::SolverIterative<double,int>*>(ksp);
  if (SS){
    info = SS->getNumIterationsTaken(this->its); CHKERRQ(info);
  }
  *iters=this->its;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverESI::SetTolerances"
int TaoLinearSolverESI::SetTolerances(double rtol, double atol, double dtol, int maxits){
  esi::SolverIterative<double,int> *SS;
  esi::ErrorCode info;
  TaoFunctionBegin;
  SS= dynamic_cast<esi::SolverIterative<double,int>*>(ksp);
  if (SS){
    info = SS->setMaxIterations( maxits ); CHKERRQ(info);
    //  info = this->ksp->setTolerance( rtol ); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverESI::View"
int TaoLinearSolverESI::View(){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverESI::SetOptions"
int TaoLinearSolverESI::SetOptions(){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

