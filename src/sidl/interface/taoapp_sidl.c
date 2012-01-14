// Implementation of wrapper class TaoSIDLApplication
// inherits from: TaoApplication

#include <iostream>
#include "sidl.hxx"
#include "Optimize_OptimizationModel.hxx"
#include "tao.h"
#include "petsc.h"
#include "src/petsctao/vector/taovec_petsc.h"
#include "src/petsctao/matrix/taomat_petsc.h"
#include "taoapp_sidl.h"
#include "LinearAlgebra_Factory.hxx"



#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::RegisterEvents"
int TaoSIDLApplication::RegisterEvents(){
  Tao_ObjectiveEval=0;
  Tao_GradientEval=0;
  Tao_HessianEval=0;
  Tao_ObjectiveAndGradientEval=0;
  Tao_Monitor = 0;
  Tao_Convergence=0;
  PetscLogEventRegister(&Tao_ObjectiveEval,"ObjEval",0);
  PetscLogEventRegister(&Tao_GradientEval,"GradEval",0);
  PetscLogEventRegister(&Tao_HessianEval,"HessEval",0);
  PetscLogEventRegister(&Tao_ObjectiveAndGradientEval,"ObjGradEval",0);
  PetscLogEventRegister(&Tao_Monitor,"Monitor",0);
  PetscLogEventRegister(&Tao_Convergence,"ConvTest",0);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::TaoSIDLApplication"
TaoSIDLApplication::TaoSIDLApplication(Optimize::OptimizationModel m)
{
  int info;
  this->taox  = this->taoxl = this->taoxu = 0;
  this->taoH = 0;
  this->setup = 0;
  this->model = m;
  this->LAFactory = 0;
  this->ksptmp=0;
  this->useh2=TAO_TRUE;
  int flag;


  // Initialize logging
  info=this->RegisterEvents();
  
  info = MPI_Initialized(&flag);
  
  if (!flag)
    return;

    
  // Call user defined initialize function
  info = this->model.initialize(); 
  this->n = this->model.getNumberVariables(); 

  Vec v;
  TaoVecPetsc *tvp;
  info = VecCreateSeq(PETSC_COMM_SELF, n, &v);
  info = TaoWrapPetscVec(v,&tvp); 
  taox = static_cast<TaoVec*>(tvp);

  Mat mat;
  TaoMatPetsc *tmp;
  info = MatCreateSeqDense(PETSC_COMM_SELF, n, n, PETSC_NULL, &mat);
  info = TaoWrapPetscMat(mat,&tmp);
  taoH = static_cast<TaoMat*>(tmp);
    

  return;
}


#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::TaoSIDLApplication-2"
TaoSIDLApplication::TaoSIDLApplication(Optimize::OptimizationModel m, LinearAlgebra::Factory factory)
{
  int info;
  int flag;
  this->taox  = this->taoxl = this->taoxu = 0;
  this->taoH = 0;
  this->setup = 0;
  this->model = m;
  this->LAFactory = factory;
  this->ksptmp=0;
  this->useh2=TAO_TRUE;

  // Initialize logging
  info=this->RegisterEvents();

  info = MPI_Initialized(&flag);
  
  if (!flag)
    return;


  // Call user defined initialize function
  info = this->model.initialize(); 
  this->n = this->model.getNumberVariables(); 


  if (LAFactory._not_nil()) {
    taox = static_cast<TaoVec*>(LAFactory.createVec(this->n));
    taoH = static_cast<TaoMat*>(LAFactory.createMat(this->n, this->n));
  } else {
    Vec v;
    TaoVecPetsc *tvp;
    info = VecCreateSeq(PETSC_COMM_SELF, n, &v);
    info = TaoWrapPetscVec(v,&tvp);
    taox = static_cast<TaoVec*>(tvp);

    Mat mat;
    TaoMatPetsc *tmp;
    info = MatCreateSeqDense(PETSC_COMM_SELF, n, n, PETSC_NULL, &mat);
    info = TaoWrapPetscMat(mat,&tmp);
    taoH = static_cast<TaoMat*>(tmp);
  }

  return;
}


#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::~TaoSIDLApplication"
TaoSIDLApplication::~TaoSIDLApplication()
{
  int info;
  info = this->model.finalize();

  info = TaoVecDestroy(taox);
  info = TaoMatDestroy(taoH);

  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::InitializeVariables"
int TaoSIDLApplication::InitializeVariables(TaoVec *xx)
{
  double *x;
  int dim,info;
  int lower = 0, upper = this->n-1, stride = 1;
  sidl::array<double> xsidl;

  info = xx->GetArray(&x,&dim); CHKERRQ(info);

  xsidl.borrow(x, 1, &lower, &upper, &stride);
  this->model.initializeVariables(xsidl);

  info = xx->RestoreArray(&x,&dim); CHKERRQ(info);

  
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::GetVariableVector"
int TaoSIDLApplication::GetVariableVector(TaoVec **xx)
{
  TaoFunctionBegin;

  if (this->taox)
    *xx = this->taox;
  else
    SETERRQ(1,"Variable TaoVec not set!");

  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::EvaluateObjectiveFunction"
int TaoSIDLApplication::EvaluateObjectiveFunction(TaoVec *xx, double *ff)
{
  // convert xx to array, call model's function
  double *x;
  int lower = 0;
  int upper = this->n-1;
  int stride = 1;
  int dim;
  sidl::array<double> xarray;
  int info;
  
  TaoFunctionBegin;
  PetscStackPush("TAO User Objective function"); 
  info = PetscLogEventBegin(Tao_ObjectiveEval,0,0,0,0);  
  info = xx->GetArray(&x,&dim); CHKERRQ(info);

  xarray.borrow(x,1,&lower,&upper,&stride);

  this->model.evaluateObjectiveFunction(xarray,*ff);

  info = xx->RestoreArray(&x,&dim); CHKERRQ(info);
  info = PetscLogEventEnd(Tao_ObjectiveEval,0,0,0,0);
  PetscStackPop;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::EvaluateGradient"
int TaoSIDLApplication::EvaluateGradient(TaoVec *xx, TaoVec *gg)
{
  double *x, *g;
  int lower = 0;
  int upper = this->n-1;
  int stride = 1;
  int dim;
  sidl::array<double> xarray;
  sidl::array<double> garray;
  int info;

  TaoFunctionBegin;
  PetscStackPush("TAO User Gradient function");
  info = PetscLogEventBegin(Tao_GradientEval,0,0,0,0);
  info = xx->GetArray(&x,&dim); CHKERRQ(info);
  info = gg->GetArray(&g,&dim); CHKERRQ(info);

  xarray.borrow(x,1,&lower,&upper,&stride);
  garray.borrow(g,1,&lower,&upper,&stride);

  this->model.evaluateGradient(xarray,garray);

  info = xx->RestoreArray(&x,&dim); CHKERRQ(info);
  info = gg->RestoreArray(&g,&dim); CHKERRQ(info);


  info = PetscLogEventEnd(Tao_GradientEval,0,0,0,0);
  PetscStackPop;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::EvaluateObjectiveAndGradient"
int TaoSIDLApplication::EvaluateObjectiveAndGradient(TaoVec *xx, double *ff, TaoVec *gg)
{
  int lower=0;
  int upper = this->n-1;
  int stride = 1;
  int dim;
  double *x, *g;
  sidl::array<double> xarray;
  sidl::array<double> garray;
  int info;
  
  TaoFunctionBegin;
  PetscStackPush("TAO User Function and Gradient function");
  info = PetscLogEventBegin(Tao_ObjectiveAndGradientEval,0,0,0,0);  

  info = xx->GetArray(&x,&dim); CHKERRQ(info);
  info = gg->GetArray(&g,&dim); CHKERRQ(info);

  xarray.borrow(x,1,&lower,&upper,&stride);
  garray.borrow(g,1,&lower,&upper,&stride);

  this->model.evaluateObjectiveAndGradient(xarray,*ff,garray);

  info = xx->RestoreArray(&x,&dim); CHKERRQ(info);
  info = gg->RestoreArray(&g,&dim); CHKERRQ(info);

  info = PetscLogEventEnd(Tao_ObjectiveAndGradientEval,0,0,0,0);
  PetscStackPop;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::EvaluateHessian"
int TaoSIDLApplication::EvaluateHessian(TaoVec *xx, TaoMat *HH)
{
  int lower[2]={0,0};
  int upper[2] = {this->n-1, this->n-1};
  int stride[2] = {1,this->n};
  int dim;
  sidl::array<double> xarray;
  sidl::array<double> harray;
  double *x_dbl, *rawh;
  TaoMatPetsc *TMPH;
  Mat H;
  MatStructure flg;
  int info=0;

  TaoFunctionBegin;
  TMPH = dynamic_cast<TaoMatPetsc*>(HH); CHKERRQ(!TMPH);
  info = TMPH->GetMatrix(&H,0,&flg); CHKERRQ(info);
  PetscStackPush("TAO User Hessian evaluation function");
  info = PetscLogEventBegin(Tao_HessianEval,0,0,0,0);
  // Get raw arrays from the tao and petsc structures
  info = xx->GetArray(&x_dbl,&dim); CHKERRQ(info);
  info = MatGetArray(H,&rawh); CHKERRQ(info);

  // wrap these raw arrays with sidl arrays
  xarray.borrow(x_dbl,1,lower,upper,stride);
  harray.borrow(rawh,2,lower,upper,stride);

  // call the user evaluation routine
  this->model.evaluateHessian(xarray,harray);
  
  // return the handles to the raw arrays
  info = MatRestoreArray(H,&rawh); CHKERRQ(info);
  info = xx->RestoreArray(&x_dbl, &dim); CHKERRQ(info);

  info = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);
  info = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);

  info = PetscLogEventEnd(Tao_HessianEval,0,0,0,0);
  PetscStackPop;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::GetHessianMatrix"
int TaoSIDLApplication::GetHessianMatrix(TaoMat **taoH)
{
  int info=0;

  TaoFunctionBegin;
  if (this->taoH)
    *taoH = (TaoMat*)this->taoH;
      
  TaoFunctionReturn(info);
}


#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::HessianSolve"
int TaoSIDLApplication::HessianSolve(TaoVec *vv, TaoVec *ww,TaoTruth *flag)
{
  double *v, *w,*x,dd=0;
  int lower = 0;
  int upper = this->n-1;
  int stride = 1;
  int dim;
  sidl::array<double> varray;
  sidl::array<double> warray;
  sidl::array<double> xarray;
  int info;
  int i;

  TaoFunctionBegin;
  PetscStackPush("TAO User Approximate Hessian function");
  info = PetscLogEventBegin(Tao_HessianEval,0,0,0,0);
  info = vv->GetArray(&v,&dim); CHKERRQ(info);
  info = ww->GetArray(&w,&dim); CHKERRQ(info);

  varray.borrow(v,1,&lower,&upper,&stride);
  warray.borrow(w,1,&lower,&upper,&stride);
  if (this->useh2==TAO_TRUE){
    TaoVec *xx0;
    info=vv->Clone(&xx0); CHKERRQ(info);
    info=TaoLMVMGetX0(tao,xx0); CHKERRQ(info);
    info = xx0->GetArray(&x,&dim); CHKERRQ(info);
    xarray.borrow(x,1,&lower,&upper,&stride);
    this->model.hessianSolve2(xarray,varray,warray);
    info = xx0->RestoreArray(&x,&dim); CHKERRQ(info);
    info=TaoVecDestroy(xx0);CHKERRQ(info);
  } else{
    this->model.hessianSolve(varray,warray);
  }
  info=vv->GetDimension(&dim);CHKERRQ(info);
  for (i=0;i<dim;i++){
    dd+=varray[i]*warray[i];
  }
  if (dd<=0){
    PetscPrintf(PETSC_COMM_WORLD,"TAOSIDLAPP: Rejecting User HessianSolve().  Inner product %4.4e <= 0\n",dd);
  }
  info = vv->RestoreArray(&v,&dim); CHKERRQ(info);
  info = ww->RestoreArray(&w,&dim); CHKERRQ(info);
  *flag=TAO_TRUE;

  info = PetscLogEventEnd(Tao_HessianEval,0,0,0,0);
  PetscStackPop;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::EvaluateVariableBounds"
int TaoSIDLApplication::EvaluateVariableBounds(TaoVec *xxll, TaoVec *xxuu)
{
  int lower[2]={0,0};
  int upper[2] = {this->n-1, this->n-1};
  int stride[2] = {1,this->n};
  int dim;
  int info;

  sidl::array<double> xl_array;
  sidl::array<double> xu_array;
  double *xl_, *xu_;

  TaoFunctionBegin;

  info = xxll->SetToConstant(TAO_NINFINITY);
  info = xxuu->SetToConstant(TAO_INFINITY);

  info = xxll->GetArray(&xl_, &dim); CHKERRQ(info);
  info = xxuu->GetArray(&xu_, &dim); CHKERRQ(info);

  // wrap these raw arrays with sidl arrays
  xl_array.borrow(xl_,1,lower,upper,stride);
  xu_array.borrow(xu_,1,lower,upper,stride);

  // Pass the sidl arrays to the model
  this->model.getVariableBounds(xl_array, xu_array); 

  info = xxll->RestoreArray(&xl_, &dim); CHKERRQ(info);
  info = xxuu->RestoreArray(&xu_, &dim); CHKERRQ(info);


  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::Monitor"
int TaoSIDLApplication::Monitor() {
  int info;
  int converged=0;
  TaoFunctionBegin;

  PetscStackPush("TAO user monitor routine");
  info = PetscLogEventBegin(Tao_Monitor,0,0,0,0);
  this->model.monitor();
  info = PetscLogEventEnd(Tao_Monitor,0,0,0,0);
  PetscStackPop;

  PetscStackPush("TAO user convergence routine");
  info = PetscLogEventBegin(Tao_Convergence,0,0,0,0);
  this->model.checkConvergence(converged);
  info = PetscLogEventEnd(Tao_Convergence,0,0,0,0);
  PetscStackPop;

  if (converged) {
    info = PetscInfo(tao,"TaoSIDLApplication::Monitor: SIDL model set convergence flag\n"); CHKERRQ(info);
    info = TaoSetTerminationReason(tao, TAO_CONVERGED_USER); CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::setTao"
int TaoSIDLApplication::setTao(TAO_SOLVER t) {
  TaoFunctionBegin;
  this->tao = t;
  TaoFunctionReturn(0);
}
  
#undef __FUNCT__
#define __FUNCT__ "TaoSIDLApplication::GetlinearSolver"
int TaoSIDLApplication::GetLinearSolver(TaoMat *H, int stype, TaoLinearSolver **SS)
{
  MPI_Comm comm2;
  int info;
  MatType mtype;
  PC pc;
  KSP ksp;
  TaoLinearSolver *S;
  PetscTruth flg1,flg2;
  Mat pm;
  MatStructure flg;
  TaoMatPetsc *TMPH;


  TaoFunctionBegin;
  TMPH = dynamic_cast<TaoMatPetsc*>(H); CHKERRQ(!TMPH);
  info = TMPH->GetMatrix(&pm,0,&flg); CHKERRQ(info);
  ksp=0;
  if (this->ksptmp){
    info = TaoWrapKSP( ksptmp, (TaoLinearSolverPetsc**)&S ); CHKERRQ(info);
    PetscObjectDereference((PetscObject)(ksptmp)); CHKERRQ(info);
    *SS=S;
    this->ksptmp=0;
    TaoFunctionReturn(0);
  } 

  if (pm){
    info = PetscObjectGetComm((PetscObject)pm,&comm2); CHKERRQ(info);
    info = MatGetType(pm,&mtype); CHKERRQ(info);
    info = PetscInfo1(pm, "mtype =%s\n",mtype); CHKERRQ(info);
    info = PetscStrncmp(mtype,MATSEQDENSE,10,&flg1); CHKERRQ(info);
    info = PetscStrncmp(mtype,MATMPIDENSE,10,&flg2); CHKERRQ(info);
  } else {
    comm2 = PETSC_COMM_WORLD;
  }
  info = KSPCreate(comm2,&ksp); CHKERRQ(info);
  info = KSPGetPC(ksp,&pc); CHKERRQ(info);
  if (flg1==PETSC_TRUE || flg2==PETSC_TRUE){
    info=PCSetType(pc,PCJACOBI); CHKERRQ(info);
  } else {
    info=PCSetType(pc,PCBJACOBI); CHKERRQ(info);
  }

  if (1 || stype==300){
    info=KSPSetType(ksp,KSPCG); CHKERRQ(info);
    info=KSPSetFromOptions(ksp); CHKERRQ(info);
  } else if (stype==310){
    info=KSPSetFromOptions(ksp); CHKERRQ(info);
    info=KSPSetType(ksp,KSPGMRES); CHKERRQ(info);
  } else if (stype==220){
    info=KSPGetPC(ksp,&pc); CHKERRQ(info);
    info=PCSetType(pc,PCNONE); CHKERRQ(info);
    info=PCSetFromOptions(pc); CHKERRQ(info);
    info=KSPSetType(ksp,KSPSTCG); CHKERRQ(info);
  } else if (stype==110){
    info=KSPSetFromOptions(ksp); CHKERRQ(info);
    info=KSPSetType(ksp,KSPGMRES); CHKERRQ(info);
  } else {
    info=KSPSetType(ksp,KSPGMRES); CHKERRQ(info);
    info=KSPSetFromOptions(ksp); CHKERRQ(info);
  }

  info = TaoWrapKSP( ksp, (TaoLinearSolverPetsc**)&S ); CHKERRQ(info);

  PetscObjectDereference((PetscObject)(ksp)); CHKERRQ(info);
  *SS=S;


  TaoFunctionReturn(0);
}







