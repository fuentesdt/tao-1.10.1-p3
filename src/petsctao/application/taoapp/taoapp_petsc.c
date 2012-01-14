#include "taoapp_petsc.h"  /*I  "tao.h"  I*/
#include "taoapp.h"
#include "src/petsctao/include/taopetsc.h"

extern int TaoVecDestroy(TaoVec*);
extern int TaoMatDestroy(TaoMat*);

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::TaoPetscApplication"
TaoPetscApplication::TaoPetscApplication(MPI_Comm mpicomm){
   this->comm=mpicomm;
   this->Setup();
   return;
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::Setup"
int TaoPetscApplication::Setup(){

  PetscFunctionBegin;

  this->taox=0;this->taoh=0;
  this->taofv=0; this->taoj=0;
  this->taofvll=0; this->taofvuu=0; this->taoAA=0;

  this->tao=0;
  this->papp=0;

  this->ksptmp = 0;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::TakeDown"
int TaoPetscApplication::TakeDown()
{
  int info;

  PetscFunctionBegin;
  info = TaoVecDestroy(taox); CHKERRQ(info);
  info = TaoMatDestroy(taoh); CHKERRQ(info);
  info = TaoVecDestroy(taofv); CHKERRQ(info);
  info = TaoMatDestroy(taoj); CHKERRQ(info);
  info = TaoVecDestroy(taofvll); CHKERRQ(info);
  info = TaoVecDestroy(taofvuu); CHKERRQ(info);
  info = TaoMatDestroy(taoAA); CHKERRQ(info);
  if (ksptmp) {
    delete ksptmp;
    ksptmp=0;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::~TaoPetscApplication"
TaoPetscApplication::~TaoPetscApplication()
{
  this->TakeDown();
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoPetscApplication::GetComm"
int TaoPetscApplication::GetComm(MPI_Comm *c) {
  PetscFunctionBegin;
  *c = comm;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::GetVariableVector"
int TaoPetscApplication::GetVariableVector(TaoVec **xx){
  int info;
  Vec V;
  PetscFunctionBegin;
  info=TaoAppGetSolutionVec(this->papp,&V); CHKERRQ(info);
  if (this->taox==0){
    info = TaoWrapPetscVec(V,&this->taox); CHKERRQ(info);
  }
  info = this->taox->SetVec(V); CHKERRQ(info);
  *xx= this->taox;
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::EvaluateObjectiveFunction"
int TaoPetscApplication::EvaluateObjectiveFunction(TaoVec *tx, double *ff)
{
  TaoVecPetsc *px = dynamic_cast <TaoVecPetsc *> (tx);

  int     info;
  PetscFunctionBegin;
  info=TaoAppComputeObjective(papp,px->GetVec(),ff);CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::EvaluateGradient"
int TaoPetscApplication::EvaluateGradient(TaoVec *tx, TaoVec *tg)
{
  TaoVecPetsc *px = dynamic_cast <TaoVecPetsc *> (tx);
  TaoVecPetsc *pg = dynamic_cast <TaoVecPetsc *> (tg);
  int     info;
  PetscFunctionBegin;
  info=TaoAppComputeGradient(papp,px->GetVec(),pg->GetVec());CHKERRQ(info);
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::EvaluateObjectiveAndGradient"
int TaoPetscApplication::EvaluateObjectiveAndGradient(TaoVec *tx, double *ff, 
						      TaoVec *tg)
{
  TaoVecPetsc *px = dynamic_cast <TaoVecPetsc *> (tx);
  TaoVecPetsc *pg = dynamic_cast <TaoVecPetsc *> (tg);
  int     info;
  PetscFunctionBegin;
  info=TaoAppComputeObjectiveAndGradient(papp,px->GetVec(),ff,pg->GetVec()); CHKERRQ(info);
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::GetHessianMatrix"
int TaoPetscApplication::GetHessianMatrix(TaoMat **HH){
  int info;
  Mat H,HP;
  MatStructure flag;

  PetscFunctionBegin;
  info=TaoAppGetHessianMat(papp,&H,&HP); CHKERRQ(info);
  if (this->taoh==0){
    info = TaoWrapPetscMat(H,&this->taoh); CHKERRQ(info);
  }
  info=this->taoh->GetMatrix(0,0,&flag);CHKERRQ(info);
  info=this->taoh->SetMatrix(H,HP,flag);CHKERRQ(info);
  *HH=this->taoh;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::EvaluateHessian"
int TaoPetscApplication::EvaluateHessian(TaoVec *tx, TaoMat *th)
{
  TaoVecPetsc *px = dynamic_cast <TaoVecPetsc *> (tx);
  TaoMatPetsc *ph = dynamic_cast <TaoMatPetsc *> (th);

  int     info;
  Mat H, HPre;
  MatStructure pflag;

  PetscFunctionBegin;
  //  info = TaoAppGetGradientVec(this,&this->FDGrad);CHKERRQ(info);
  info = ph->GetMatrix(&H,&HPre,&pflag); CHKERRQ(info);
  info = TaoAppComputeHessian(this->papp, px->GetVec(), &H, &HPre, &pflag); CHKERRQ(info);
  info = ph->SetMatrix(H,HPre,pflag); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::HessianSolve"
int TaoPetscApplication::HessianSolve(TaoVec *tv, TaoVec *tw, TaoTruth *success)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *pw = dynamic_cast <TaoVecPetsc *> (tw);

  int     info;
  PetscTruth success2;
  PetscFunctionBegin;
  info=TaoAppHessianSolve(papp,pv->GetVec(),pw->GetVec(),&success2);CHKERRQ(info);
  if (success2==PETSC_TRUE){*success=TAO_TRUE;} else { *success=TAO_FALSE;}
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::EvaluateVariableBounds"
int TaoPetscApplication::EvaluateVariableBounds(TaoVec *tl, TaoVec *tu)
{
  TaoVecPetsc *pl = dynamic_cast <TaoVecPetsc *> (tl);
  TaoVecPetsc *pu = dynamic_cast <TaoVecPetsc *> (tu);
  int info;
  PetscFunctionBegin;
  info=TaoAppComputeVariableBounds(this->papp,pl->GetVec(),pu->GetVec()); CHKERRQ(info);
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::GetConstraintVector"
int TaoPetscApplication::GetConstraintVector(TaoVec **rr){
  int info;
  Vec R;
  PetscFunctionBegin;
  info=TaoAppGetFunctionVec(this->papp,&R); CHKERRQ(info);
  if (this->taofv==0){
    info = TaoWrapPetscVec(R,&this->taofv); CHKERRQ(info);
  }
  info = this->taofv->SetVec(R); CHKERRQ(info);
  *rr= this->taofv;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::EvaluateConstraints"
int TaoPetscApplication::EvaluateConstraints(TaoVec *tx, TaoVec *tr)
{
  TaoVecPetsc *px = dynamic_cast <TaoVecPetsc *> (tx);
  TaoVecPetsc *pr = dynamic_cast <TaoVecPetsc *> (tr);
  int     info;

  PetscFunctionBegin;
  info=TaoAppComputeFunction(papp,px->GetVec(),pr->GetVec());CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::GetJacobian"
int TaoPetscApplication::GetJacobianMatrix(TaoMat **JJ){
  int info;
  Mat J,JP;
  MatStructure flag;

  PetscFunctionBegin;
  info=TaoAppGetJacobianMat(this->papp,&J,&JP); CHKERRQ(info);
  if (this->taoj==0){
    info = TaoWrapPetscMat(J,&this->taoj); CHKERRQ(info);
  }
  info=this->taoj->GetMatrix(0,0,&flag);CHKERRQ(info);
  info=this->taoj->SetMatrix(J,JP,flag);CHKERRQ(info);
  *JJ=this->taoj;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::EvaluateJacobian"
int TaoPetscApplication::EvaluateJacobian(TaoVec *tx, TaoMat *tj)
{
  TaoVecPetsc *px = dynamic_cast <TaoVecPetsc *> (tx);
  TaoMatPetsc *pj = dynamic_cast <TaoMatPetsc *> (tj);

  int     info;
  Mat J,JPre;
  MatStructure pflag;

  PetscFunctionBegin;
  info = pj->GetMatrix(&J,&JPre,&pflag);CHKERRQ(info);
  info = TaoAppComputeJacobian(papp,px->GetVec(),&J,&JPre,&pflag);CHKERRQ(info);
  info = pj->SetMatrix(J,JPre,pflag);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::GetInequalityConstraints"
int TaoPetscApplication::GetInequalityConstraints(TaoVec **LL, TaoMat **JJ, TaoVec**UU){
  PetscFunctionBegin;
  *LL=this->taofvll;
  *JJ=this->taoAA;
  *UU=this->taofvuu;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoPetscApplication::GetLinearSolver"
int TaoPetscApplication::GetLinearSolver(TaoMat *tH, int stype, TaoLinearSolver **SS)
{
  TaoMatPetsc *pH = dynamic_cast <TaoMatPetsc *> (tH);

  Mat pm, pmpre;
  MatStructure pflag;

  MPI_Comm comm2;
  int info;
  const MatType mtype;
  PC pc;
  KSP newksp;
  TaoLinearSolverPetsc *S;
  PetscTruth flg1,flg2;

  PetscFunctionBegin;

  info = pH->GetMatrix(&pm,&pmpre,&pflag); CHKERRQ(info);

  info = TaoAppGetKSP(papp,&newksp); CHKERRQ(info);
  if (this->ksptmp==0){
    info = TaoWrapKSP( newksp, &S ); CHKERRQ(info);
    this->ksptmp=S;
  }

  info=this->ksptmp->SetKSP(newksp);CHKERRQ(info);
  *SS=this->ksptmp;

  if (pm) {
    info = PetscObjectGetComm((PetscObject)pm,&comm2); CHKERRQ(info);
    info = MatGetType(pm,&mtype); CHKERRQ(info);
    info = PetscStrncmp(mtype,MATSEQDENSE,10,&flg1); CHKERRQ(info);
    info = PetscStrncmp(mtype,MATMPIDENSE,10,&flg2); CHKERRQ(info);
  } 
  else {
    comm2 = this->comm;
  }

  info = KSPGetPC(newksp,&pc); CHKERRQ(info);
  info=KSPSetFromOptions(newksp); CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::InitializeVariables"
int TaoPetscApplication::InitializeVariables(TaoVec *xx){
  PetscFunctionBegin; 
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::Monitor"
int TaoPetscApplication::Monitor(){
  int info;
  PetscFunctionBegin;
  info=TaoAppMonitor(this->papp); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplication::Monitor2"
int TaoPetscApplication::Monitor2(TaoVec *xx, TaoVec *tg, TaoVec *dx, TaoTruth *flag){
  int     info;
  PetscTruth flag2;

  PetscFunctionBegin;
  *flag = TAO_FALSE;;
  if (tg) {
    TaoVecPetsc *pg = dynamic_cast <TaoVecPetsc *> (tg);
    info=TaoAppCheckConvergence(papp, pg->GetVec(), &flag2); CHKERRQ(info);
    if (flag2==PETSC_TRUE) *flag=TAO_TRUE;
  } 
  PetscFunctionReturn(0);
}


