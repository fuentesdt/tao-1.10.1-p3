#include "tao_solver.h"
#include "lm.h"

int Monitor_NLSQ(TAO_SOLVER tao, void *ctx);
extern int TaoSetMaxGPIts(TAO_SOLVER tao, int its);


//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoSetup_NLSQ"
int TaoSetup_NLSQ(TAO_SOLVER tao, void *ctx)
{
  TAO_NLSQ *nlsq = (TAO_NLSQ*)ctx;
  int info;
  char type[256];
  TaoVec *XX,*RR,*LL,*UU;
  TaoMat *JJ;
  TaoLinearSolver *ksp;
  int dim;  // Dimension of XX
  int factor=3;  // multiple of dim to use as default ksp_max_its
  TaoTruth flg;

  TaoFunctionBegin;

  nlsq->lsq = new TaoLeastSquaresApplication();

  
  info = TaoCheckConstraints(tao); CHKERRQ(info);

  //             TaoVec GG
  //    Getting Handles to existing objects nlsq->JJ, nlsq->RR 

  // one dimension is [x]=n, one dimension is [RR=c(x)]=m
  //   so -- get these vectors for cloning

  info = TaoGetSolution(tao,&XX); CHKERRQ(info);   // R^n
  
  // The matrix RR will point to is already created by the user (should check first, though)
  info = TaoGetConstraints(tao,&RR); CHKERRQ(info);
  info = TaoGetJacobian(tao,&JJ); CHKERRQ(info);

  info = TaoComputeJacobian(tao,XX,JJ);CHKERRQ(info);
  info = nlsq->lsq->SetupCG(XX,RR,JJ); CHKERRQ(info);
  nlsq->lsq->outertao=tao;
  // Create the new trust region solver
  //set up bounds
  info = TaoGetVariableBounds(tao,&LL, &UU); CHKERRQ(info);
  if (LL && UU){
    info = TaoCheckBounds(tao); CHKERRQ(info);
    info = TaoCreateFull("tao_tron","nlsq_",tao->comm,&nlsq->innertao);
    CHKERRQ(info);
    info = TaoSetMaxGPIts(nlsq->innertao,1); CHKERRQ(info);
  } else {
    info = TaoCreateFull("tao_ntr","nlsq_",tao->comm,&nlsq->innertao);
    CHKERRQ(info);
  }
  info = TaoOptionString("-lsqr_tao_method","The underlying minimization solver type ","","",type,256, &flg);CHKERRQ(info);
  if (flg) {
    info = TaoSetMethod(nlsq->innertao,(TaoMethod) type);CHKERRQ(info);
  }

  //  info = TaoGetLinearSolver(tao,&ksp); CHKERRQ(info);
  //  info = TaoSetLinearSolver(nlsq->innertao,ksp); CHKERRQ(info);

  // Set up default criteria
  info= TaoSetApplication(nlsq->innertao, nlsq->lsq); CHKERRQ(info);
  info = XX->GetDimension(&dim); CHKERRQ(info);
  info = TaoGetLinearSolver(nlsq->innertao,&ksp); CHKERRQ(info);
  if (ksp){
    info = ksp->SetTolerances(1.0e-4,1.0e-30,1.0e30,dim*factor);CHKERRQ(info);
    info = TaoSetLinearSolver(tao,ksp); CHKERRQ(info);
  }

  if (tao->vec_grad==0){ tao->vec_grad=nlsq->lsq->GG;}

  tao->hessian=nlsq->lsq->HH;
  tao->vec_sol_update=nlsq->innertao->vec_sol_update;
  
  TaoFunctionReturn(0);
}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoDestroy_NLSQ"
int TaoDestroy_NLSQ(TAO_SOLVER tao, void *solver)
{
  TAO_NLSQ *nlsq = (TAO_NLSQ*)solver;
  int info;

  TaoFunctionBegin;

  if (tao->setupcalled) {

    // Destroy the trust-region solver
    info = TaoDestroy(nlsq->innertao); CHKERRQ(info);
    info = TaoApplicationDestroy(nlsq->lsq); CHKERRQ(info);
    info = TaoDestroyLinearSolver(tao); CHKERRQ(info);

    // Destroy any vectors or matrices that were allocated
    info = TaoMatDestroy(nlsq->lsq->HH); CHKERRQ(info);

    if (nlsq->lsq->JJ){
      info = TaoMatDestroy(nlsq->lsq->JJ); CHKERRQ(info);
    }
    info = TaoVecDestroy(nlsq->lsq->DD1); CHKERRQ(info);
    info = TaoVecDestroy(nlsq->lsq->DD2); CHKERRQ(info);
    info = TaoVecDestroy(nlsq->lsq->GG); CHKERRQ(info);
  
  }
  // Finally, destroy the lm context

  info = TaoFree(nlsq); CHKERRQ(info);
  

  TaoFunctionReturn(0);
}


//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoSolve_NLSQ"
  int TaoSolve_NLSQ(TAO_SOLVER tao, void *solver)
{
  TAO_NLSQ *nlsq = (TAO_NLSQ*)solver;
  int info;
  TaoVec *XX,*DX;

  TaoFunctionBegin;

  info = TaoGetSolution(tao,&XX); CHKERRQ(info);
  info = TaoGetSolutionUpdate(nlsq->innertao,&DX); CHKERRQ(info);
  tao->vec_sol_update=DX;
  info = TaoSetTrustRegionRadius(nlsq->innertao,tao->trust0); CHKERRQ(info);
  info = TaoClearMonitor(nlsq->innertao); CHKERRQ(info);
  // Do not use innertao's convergence test
  info = TaoSetConvergenceTest(nlsq->innertao,TAO_NULL,TAO_NULL); CHKERRQ(info);
  info = TaoSetMonitor(nlsq->innertao,Monitor_NLSQ,(void*)tao,TAO_NULL); CHKERRQ(info);
  info = TaoSolve(nlsq->innertao); CHKERRQ(info);

  TaoFunctionReturn(0);
}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoSetOptions_NLSQ"
int TaoSetOptions_NLSQ(TAO_SOLVER tao, void *solver)
{
  int info;
  char type[20];
  TAO_NLSQ *nlsq = (TAO_NLSQ*)solver;
  TAO_SOLVER innertao=nlsq->innertao;
  TaoTruth flg;

  TaoFunctionBegin;
  info = TaoOptionString("-lsqr_tao_method","The underlying minimization solver type ","","tao_ntr or tao_tron",type,20, &flg);CHKERRQ(info);
  if (innertao->setfromoptions) {
    info = (*innertao->setfromoptions)(innertao,innertao->data); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}


//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoView_NLSQ"
  int TaoView_NLSQ(TAO_SOLVER tao, void *solver)
{
  TAO_NLSQ *nlsq = (TAO_NLSQ*)solver;
  int info;
  TaoFunctionBegin;
  info=TaoPrintStatement(tao,"Minimization Solver --\n");CHKERRQ(info);
  info=TaoPrintStatement(tao,"  Ignore Convergence Tolerances: \n");CHKERRQ(info);
  info=TaoView(nlsq->innertao);CHKERRQ(info);
  info=TaoPrintStatement(tao,"-- \n");CHKERRQ(info);
  TaoFunctionReturn(0);
}

//========================================
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "TaoCreate_NLSQ"
int TaoCreate_NLSQ(TAO_SOLVER tao)
{
  TAO_NLSQ *nlsq;
  int info;
  double fatol = 1e-8, frtol =1e-12, catol = 0.0, crtol = 0.0;  

  TaoFunctionBegin;  

  // Allocate the tao_lm context
  info = TaoNew(TAO_NLSQ, &nlsq); CHKERRQ(info);
  
  info = TaoSetSolver(tao, TaoSetup_NLSQ, TaoSetOptions_NLSQ, TaoSolve_NLSQ, 
		      TaoView_NLSQ, TaoDestroy_NLSQ,(void*)nlsq); 
  CHKERRQ(info);
  
  // Set the nlsqtao field in the context to this solver (we need this solver to calculate funcgrad, jacobian)
  //  nlsq->nlsqtao = tao;  
  info = TaoSetMaximumIterates(tao,10000); CHKERRQ(info);
  info = TaoSetTolerances(tao, fatol, frtol, catol, crtol); CHKERRQ(info);  
  info = TaoSetFunctionLowerBound(tao,0); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::TaoLeastSquaresApplication()"
TaoLeastSquaresApplication::TaoLeastSquaresApplication(){
  this->JJ=0;
  this->JJac=0;
  this->HH=0;
  this->DD1=0;
  this->DD2=0;
  this->RR=0;
  this->GG=0;
  return;
}


#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::~TaoLeastSquaresApplication()"
TaoLeastSquaresApplication::~TaoLeastSquaresApplication(){
  return;
}


#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::SetupCG"
int TaoLeastSquaresApplication::SetupCG(TaoVec*XX, TaoVec *RR, TaoMat*JJac){
  int info;
  TaoFunctionBegin;

  this->RR=RR; this->JJac=JJac; this->JJ=0;

  info = XX->Clone(&this->DD2); CHKERRQ(info); // R^m
  info = XX->Clone(&this->GG); CHKERRQ(info); // R^m  (will be JJ^T * RR) 
  info = RR->Clone(&this->DD1); CHKERRQ(info); // R^n
  
  // -- ADA Matrices manipulate Vec's and Mat's (not TaoVec's, TaoMat's), so need handles
  info = this->DD1->SetToConstant(1.0); CHKERRQ(info);
  info = this->DD2->SetToZero(); CHKERRQ(info);

  info = this->JJac->Clone(&this->JJ); CHKERRQ(info);
  info = this->JJ->CreateATDAMatrix(this->DD1,this->DD2, &this->HH); CHKERRQ(info);

  // Set nlsq->DD1 to all ones, nlsq->DD2 to all zeros (These are used to create H)
  
  //Create JJ (R^{mxn}) and HH matrices (R^{nxn}). (HH = JJ^T * diag(DD1) * J + DD2) 
  // -- JJ already created, just need a handle to it

  TaoFunctionReturn(0);
}


//========================================
#undef __FUNCT__
#define __FUNCT__ "FormFunction_NLSQ"
int  TaoLeastSquaresApplication::EvaluateObjectiveFunction(TaoVec *XX, double *fPtr)
{
  /* The function value of the trust region problem is ||r||_2^2 / 2,
     where r (already in RR) is the function of the residuals [ c(x) ]*/
  int info;

  TaoFunctionBegin;

  // Get the residual vector at this xx
  info = TaoComputeConstraints(this->outertao,XX,this->RR); CHKERRQ(info);
  // compute f = ||r||_2^2 / 2
  info = this->RR->Norm2squared(fPtr); CHKERRQ(info);
  *fPtr/=2.0;

  TaoFunctionReturn(0);
}


//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::EvaluateObjectiveAndGradient"
int TaoLeastSquaresApplication::EvaluateObjectiveAndGradient(TaoVec *XX, double *fPtr, TaoVec *GGVec)
{
  /* The function value of the trust region problem is ||r||_2^2 / 2,
     where r (already in RR) is the function of the residuals [ c(x) ]*/
  int info;

  TaoFunctionBegin;

  // Get the residual vector at this xx
  info =   TaoComputeConstraints(this->outertao,XX,this->RR); CHKERRQ(info);
  // compute f = ||r||_2^2 / 2
  info = this->RR->Norm2squared(fPtr); CHKERRQ(info);
  *fPtr/=2.0;

  /* The gradient is J^T*r */
  // Compute the jacobian
  info = TaoComputeJacobian(this->outertao,XX,this->JJac);CHKERRQ(info);
  // Calculate the gradient (G = J^t * R)
  info = this->JJac->MultiplyTranspose(this->RR,GGVec); CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::EvaluateGradient"
int TaoLeastSquaresApplication::EvaluateGradient(TaoVec *XX, TaoVec *GG)
{
  /* The function value of the trust region problem is ||r||_2^2 / 2,
     where r (already in RR) is the function of the residuals [ c(x) ]*/
  int info;
  double fPtr;

  TaoFunctionBegin;

  // Get the residual vector at this xx
  info =   TaoComputeConstraints(this->outertao,XX,this->RR); CHKERRQ(info);
  // compute f = ||r||_2^2 / 2
  info = this->RR->Norm2squared(&fPtr); CHKERRQ(info);
  fPtr/=2.0;

  /* The gradient is J^T*r */
  // Compute the jacobian
  info = TaoComputeJacobian(this->outertao,XX,this->JJac);CHKERRQ(info);
  // Calculate the gradient (G = J^t * R)
  info = this->JJac->MultiplyTranspose(this->RR,GG); CHKERRQ(info);

  TaoFunctionReturn(0);
}


//========================================
#undef __FUNCT__
#define __FUNCT__ "ComputeHessian_NLSQ"
int TaoLeastSquaresApplication::EvaluateHessian(TaoVec *XX, TaoMat *HH)
{
  /*  The hessian for the trust region problem is 2*(J^T)*J, where J is the jacobian of the constraint function */
  int info;
  TaoFunctionBegin;

  info = this->JJ->CopyFrom(this->JJac); CHKERRQ(info);
  info = this->DD2->SetToConstant(0.0); CHKERRQ(info);
  info = this->DD1->SetToConstant(1.0); CHKERRQ(info);

  // Note that since HH contains a pointer to H, which has a pointe to J, changing JJ
  //         also changes HH.  Also, JJ was updated in TaoLeastSquaresApplication::EvaluateFunctionGradient
  // This means that the hessian is already formed, so we do nothing here.
 
  TaoFunctionReturn(0);
}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::GetVariableVector"
int TaoLeastSquaresApplication::GetVariableVector(TaoVec **XX)
{
  int info;
  TaoFunctionBegin;
  info = TaoGetSolution(this->outertao,XX);CHKERRQ(info);
  TaoFunctionReturn(0);
}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::GetVariableVector"
int TaoLeastSquaresApplication::GetGradientVector(TaoVec **GVec)
{
  TaoFunctionBegin;
  if (GVec){
    *GVec=this->GG;
  }
  TaoFunctionReturn(0);
}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::GetVariableVector"
int TaoLeastSquaresApplication::GetHessianMatrix(TaoMat **HHes){

  TaoFunctionBegin;
  if (HHes){
    *HHes=this->HH;
  }
  TaoFunctionReturn(0);
}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::intializeVariables"
int TaoLeastSquaresApplication::InitializeVariables(TaoVec *XX0)
{
  int info;
  TaoFunctionBegin;
  info = this->outertao->taoappl->InitializeVariables(XX0); CHKERRQ(info);
  TaoFunctionReturn(0);
}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoLeastSquaresApplication::GetVariableBounds"
int TaoLeastSquaresApplication::GetVariableBounds(TaoVec **XXLL, TaoVec **XXUU)
{
  int info;
  TaoFunctionBegin;
  info = TaoGetVariableBounds(this->outertao,XXLL,XXUU);CHKERRQ(info);
  TaoFunctionReturn(0);
}


//========================================
#undef __FUNCT__
#define __FUNCT__ "LM_Monitor"
/* This function receives monitor information from the trust region solver and transfers the 
   data to any monitors listening to the levenberg-marquardt solver. */
int Monitor_NLSQ(TAO_SOLVER tao, void *ctx)
{
  //  TAO_NLSQ *nlsq = (TAO_NLSQ*)ctx;
  TAO_SOLVER outertao=(TAO_SOLVER)ctx;
  int iterate;
  double f,fnorm,cnorm,step;
  TaoTerminateReason reason;
  int info;


  TaoFunctionBegin;
  info = TaoGetIterationData(tao,&iterate,&f,&fnorm,&cnorm,&step,&reason); CHKERRQ(info);
  info = TaoMonitor(outertao,iterate,2*f,2*fnorm,cnorm,step,&reason); CHKERRQ(info); 
  tao->reason = reason;
  TaoFunctionReturn(0);
}
