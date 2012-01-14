#include "src/bound/impls/kt/kt.h"
#include "src/tao_impl.h"

static int TaoSetDown_KT(TAO_SOLVER, void *);
static int TaoMonitor_KT(TAO_SOLVER my_tao, void *solver);

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_KT"
int TaoSetUp_KT(TAO_SOLVER tao,void *solver)
{
  TAO_KT *kt = (TAO_KT *)solver;
  TaoVec *xx;
  TaoTruth flg;
  int info;

  TaoFunctionBegin;
  info = TaoGetSolution(tao,&xx); CHKERRQ(info);

  // Create new complementarity solver
  info = TaoOptionString("-tao_kt_method",0,0,kt->comp_method,
		         kt->comp_method,256,&flg);
  CHKERRQ(info);
  if (kt->csolver==0){
      info = TaoCreateFull(kt->comp_method,"t",((PetscObject)tao)->comm, &kt->csolver); 
    CHKERRQ(info);
  }
  kt->ktapp = new TaoKTApplication(tao,kt->csolver);
  kt->ktapp->SetItUp1(); CHKERRQ(info);
  kt->ktapp->SetItUp2(); CHKERRQ(info);
  kt->setupcalled=1;
  
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_KT"
static int TaoSetDown_KT(TAO_SOLVER tao, void *solver)
{
  TAO_KT *kt = (TAO_KT *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoDestroy(kt->csolver); CHKERRQ(info); 
  info = TaoDestroyApplication(kt->ktapp); CHKERRQ(info); 
  kt->setupcalled=0;
  kt->csolver=0;
  kt->ktapp=0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMonitor_KT"
static int TaoMonitor_KT(TAO_SOLVER complementaritytao, void *solver)
{
  TaoKTApplication* ktapp= (TaoKTApplication*)solver;
  TaoTerminateReason reason;
  double f, fnorm, cnorm, xdiff;
  TaoInt iterate;
  int info;

  TaoFunctionBegin;
  info = TaoGetSolutionStatus(complementaritytao, &iterate, &f, &fnorm, 
			     &cnorm, &xdiff, &reason); CHKERRQ(info);
  f=ktapp->func;
  info = TaoMonitor(ktapp->orig,iterate,f,fnorm,cnorm,xdiff,&reason); CHKERRQ(info);
  info = TaoSetTerminationReason(complementaritytao,reason); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "TaoSolve_KT"
static int TaoSolve_KT(TAO_SOLVER tao, void *solver)
{
  TAO_KT *kt = (TAO_KT *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoSolve(kt->csolver);  CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_KT"
static int TaoSetOptions_KT(TAO_SOLVER tao, void*solver)
{
  TAO_KT *kt = (TAO_KT *)solver;
  TAO_SOLVER csolver=kt->csolver;
  TaoTruth flg;
  int info;

  TaoFunctionBegin;
  info = TaoOptionString("-tao_kt_method",
                         "Set method for solving kt conditions",
		         "TaoKTSetMethod",kt->comp_method,
		         kt->comp_method,256,&flg);
  CHKERRQ(info);

  if (csolver==0){
      info = TaoCreateFull(kt->comp_method,"t",((PetscObject)tao)->comm, &kt->csolver); 
    CHKERRQ(info);
  }
  csolver=kt->csolver;
  if (csolver->setfromoptions) {
    info = (*csolver->setfromoptions)(csolver,csolver->data); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_KT"
static int TaoView_KT(TAO_SOLVER tao,void* solver)
{
  TAO_KT   *kt = (TAO_KT *)solver;
  int      info;

  TaoFunctionBegin;
  info = TaoPrintString(tao,"  kt method=%s\n", kt->comp_method); 
         CHKERRQ(info);
  TaoFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_KT"
int TaoCreate_KT(TAO_SOLVER tao)
{
  TAO_KT *kt;
  int    info;

  TaoFunctionBegin;

  info = TaoNew(TAO_KT,&kt); CHKERRQ(info);

  kt->csolver = 0;
  kt->setupcalled=0;
  info = TaoStrcpy(kt->comp_method, "tao_ssils"); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_KT,(void*)kt); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_KT,TaoSetDown_KT); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetOptions_KT); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_KT); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,2000); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,4000); CHKERRQ(info);

  info = TaoSetTolerances(tao,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,1.0e-12,0.0,0.0); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::TaoKTApplication"
TaoKTApplication::TaoKTApplication(TAO_SOLVER outertao, TAO_SOLVER innertao) 
{
  this->orig=outertao;
  this->csolver=innertao;
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::SetItUp1"
int TaoKTApplication::SetItUp1() 
{
  int info;
  TaoFunctionBegin;
  info = TaoSetConvergenceTest(this->csolver, TAO_NULL, TAO_NULL); CHKERRQ(info);
  info = TaoSetTolerances(this->csolver,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(this->csolver,0.0,0.0,0.0); CHKERRQ(info);
  info = TaoSetMaximumIterates(this->csolver,0); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(this->csolver,0); CHKERRQ(info);
  info = TaoSetFunctionLowerBound(this->csolver,0); CHKERRQ(info);
  info = TaoSetMonitor(this->csolver, TaoMonitor_KT, (void*)this); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::SetItUp2"
int TaoKTApplication::SetItUp2() 
{
  int info;
  TaoVec *x1,*x2;
  TaoFunctionBegin;
  info = TaoSetApplication(this->csolver,this); CHKERRQ(info);
  info = TaoGetVariableBounds(this->csolver,&x1,&x2); CHKERRQ(info);
  info = TaoSetVariableBounds(this->orig,x1,x2); CHKERRQ(info);
  info = TaoGetStepDirectionVector(this->csolver,&x1); CHKERRQ(info);
  info = TaoSetStepDirectionVector(this->orig,x1); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::EvaluateConstraints"
int TaoKTApplication::EvaluateConstraints(TaoVec *x, TaoVec *f) 
{
  int info;
  TaoFunctionBegin;
  info = TaoComputeFunctionGradient(this->orig, x, &this->func, f); 
         CHKERRQ(info);
  info = TaoSetLagrangianGradientVector(this->orig,f); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::EvaluateJacobian"
int TaoKTApplication::EvaluateJacobian(TaoVec *x, TaoMat *J) 
{
  int info;
  TaoFunctionBegin;
  info = TaoComputeHessian(this->orig, x, J); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::InitializeVariables"
int TaoKTApplication::InitializeVariables(TaoVec *x) 
{
  int info;
  TaoApplication* theapp;
  TaoFunctionBegin;
  info = TaoGetApplication(this->orig, &theapp); CHKERRQ(info);
  info = theapp->InitializeVariables(x); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::GetVariableVector"
int TaoKTApplication::GetVariableVector(TaoVec **x) 
{
  int info;
  TaoFunctionBegin;
  info = TaoGetSolution(this->orig, x); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::EvaluateVariableBounds"
int TaoKTApplication::EvaluateVariableBounds(TaoVec *l, TaoVec *u) 
{
  int info;
  TaoFunctionBegin;
  info = TaoEvaluateVariableBounds(this->orig, l, u); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::GetJacobianMatrix"
int TaoKTApplication::GetJacobianMatrix(TaoMat **J) 
{
  int info;
  TaoFunctionBegin;
  info = TaoGetHessian(this->orig, J); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoKTApplication::GetLinearSolver"
int TaoKTApplication::GetLinearSolver(TaoMat *H, TaoInt stype, TaoLinearSolver **tksp) {
  int info;
  TaoFunctionBegin;
  info = TaoCreateLinearSolver(this->orig, H, stype,tksp); CHKERRQ(info);
  TaoFunctionReturn(0);
}
