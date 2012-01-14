#include "src/complementarity/impls/ssls/ssls.h"

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_SSLS"
int TaoSetUp_SSLS(TAO_SOLVER tao, void *solver)
{
  TAO_SSLS *ssls = (TAO_SSLS *)solver;
  TaoVec *x;
  TaoMat *J;
  TaoInt n;
  int info;

  TaoFunctionBegin;

  info = TaoGetSolution(tao, &x); CHKERRQ(info);
  info = TaoGetJacobian(tao, &J); CHKERRQ(info);

  info = x->GetDimension(&n); CHKERRQ(info);

  info = x->Clone(&(ssls->f)); CHKERRQ(info);
  info = x->Clone(&(ssls->g)); CHKERRQ(info);
  info = x->Clone(&(ssls->ff)); CHKERRQ(info);
  info = x->Clone(&(ssls->dpsi)); CHKERRQ(info);
  info = x->Clone(&(ssls->d)); CHKERRQ(info);
  info = x->Clone(&(ssls->w)); CHKERRQ(info);
  info = x->Clone(&(ssls->da)); CHKERRQ(info);
  info = x->Clone(&(ssls->db)); CHKERRQ(info);
  info = x->Clone(&(ssls->t1)); CHKERRQ(info);
  info = x->Clone(&(ssls->t2)); CHKERRQ(info);
  info = x->Clone(&(ssls->xl)); CHKERRQ(info);
  info = x->Clone(&(ssls->xu)); CHKERRQ(info);

  info = TaoSetVariableBounds(tao,ssls->xl,ssls->xu); CHKERRQ(info);
  info = TaoCheckBounds(tao); CHKERRQ(info);

  info = TaoSetStepDirectionVector(tao,ssls->d);CHKERRQ(info);
  info = TaoSetLagrangianGradientVector(tao,ssls->g);CHKERRQ(info);

  info = TaoCreateLinearSolver(tao,J,210,0); CHKERRQ(info);

  info = TaoLineSearchSetUp(tao);CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_SSLS"
int TaoSetDown_SSLS(TAO_SOLVER tao, void *solver)
{
  TAO_SSLS *ssls = (TAO_SSLS *)solver;
  int info;

  TaoFunctionBegin;

  info=TaoVecDestroy(ssls->f); CHKERRQ(info);
  info=TaoVecDestroy(ssls->g); CHKERRQ(info);
  info=TaoVecDestroy(ssls->ff); CHKERRQ(info);
  info=TaoVecDestroy(ssls->dpsi); CHKERRQ(info);
  info=TaoVecDestroy(ssls->d); CHKERRQ(info);
  info=TaoVecDestroy(ssls->w); CHKERRQ(info);
  info=TaoVecDestroy(ssls->da); CHKERRQ(info);
  info=TaoVecDestroy(ssls->db); CHKERRQ(info);
  info=TaoVecDestroy(ssls->t1); CHKERRQ(info);
  info=TaoVecDestroy(ssls->t2); CHKERRQ(info);
  info=TaoVecDestroy(ssls->xl); CHKERRQ(info);
  info=TaoVecDestroy(ssls->xu); CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_SSLS"
int TaoSetOptions_SSLS(TAO_SOLVER tao, void *solver)
{
  TAO_SSLS *ssls = (TAO_SSLS *)solver;
  int info;
  TaoTruth flg;

  TaoFunctionBegin;
  info = TaoOptionsHead("Semismooth method with a linesearch for "
  		        "complementarity problems"); CHKERRQ(info);
  info = TaoOptionDouble("-ssls_delta", "descent test fraction", "",
                         ssls->delta, &(ssls->delta), &flg);CHKERRQ(info);
  info = TaoOptionDouble("-ssls_rho", "descent test power", "",
                         ssls->rho, &(ssls->rho), &flg);CHKERRQ(info);
  info = TaoLineSearchSetFromOptions(tao);CHKERRQ(info);
  info = TaoOptionsTail(); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_SSLS"
int TaoView_SSLS(TAO_SOLVER tao, void *solver)
{
  int info;

  TaoFunctionBegin;
  info = TaoLineSearchView(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "Tao_SSLS_Function"
int Tao_SSLS_Function(TAO_SOLVER tao, TaoVec *X, double *fcn, void *solver) 
{
  TAO_SSLS *ssls = (TAO_SSLS *)solver;
  TaoVec *f, *l, *u, *ff;
  int info;

  TaoFunctionBegin;
  info = TaoGetVariableBounds(tao, &l, &u); CHKERRQ(info);
  f  = ssls->f;
  ff = ssls->ff;
  
  info = TaoComputeConstraints(tao, X, f); CHKERRQ(info);
  info = ff->Fischer(X, f, l, u); CHKERRQ(info);
  info = ff->Norm2(&ssls->merit); CHKERRQ(info);
  *fcn = 0.5*ssls->merit*ssls->merit;
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "Tao_SSLS_FunctionGradient"
int Tao_SSLS_FunctionGradient(TAO_SOLVER tao, TaoVec *X, double *fcn, 
                              TaoVec *G, void *solver)
{
  TAO_SSLS *ssls = (TAO_SSLS *)solver;
  TaoVec *f, *l, *u, *ff, *da, *db, *t1, *t2;
  TaoMat *J;
  int info;

  TaoFunctionBegin;
  info = TaoGetVariableBounds(tao, &l, &u); CHKERRQ(info);
  info = TaoGetJacobian(tao, &J); CHKERRQ(info);
  f  = ssls->f;
  ff = ssls->ff;
  da = ssls->da;
  db = ssls->db;
  t1 = ssls->t1;
  t2 = ssls->t2;

  info = TaoComputeConstraints(tao, X, f); CHKERRQ(info);
  info = ff->Fischer(X, f, l, u); CHKERRQ(info);
  info = ff->Norm2(&ssls->merit); CHKERRQ(info);
  *fcn = 0.5*ssls->merit*ssls->merit;

  info = TaoComputeJacobian(tao, X, J); CHKERRQ(info);
  info = J->D_Fischer(X, f, l, u, t1, t2, da, db); CHKERRQ(info);
  info = J->RowScale(db); CHKERRQ(info);
  info = J->AddDiagonal(da); CHKERRQ(info);
  info = J->MultiplyTranspose(ff, G); CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "Tao_ASLS_FunctionGradient"
int Tao_ASLS_FunctionGradient(TAO_SOLVER tao, TaoVec *X, double *fcn, 
                              TaoVec *G, void *solver)
{
  TAO_SSLS *asls = (TAO_SSLS *)solver;
  TaoVec *f, *l, *u, *ff, *da, *db, *t1, *t2;
  TaoMat *J;
  int info;

  TaoFunctionBegin;
  info = TaoGetVariableBounds(tao, &l, &u); CHKERRQ(info);
  info = TaoGetJacobian(tao, &J); CHKERRQ(info);
  f  = asls->f;
  ff = asls->ff;
  da = asls->da;
  db = asls->db;
  t1 = asls->t1;
  t2 = asls->t2;

  info = TaoComputeConstraints(tao, X, f); CHKERRQ(info);
  info = ff->Fischer(X, f, l, u); CHKERRQ(info);
  info = ff->Norm2(&asls->merit); CHKERRQ(info);
  *fcn = 0.5*asls->merit*asls->merit;

  info = TaoComputeJacobian(tao, X, J); CHKERRQ(info);
  info = J->D_Fischer(X, f, l, u, t1, t2, da, db); CHKERRQ(info);
  info = t1->PointwiseMultiply(ff, db);
  info = J->MultiplyTranspose(t1, G); CHKERRQ(info);
  info = t1->PointwiseMultiply(ff, da); CHKERRQ(info);
  info = G->Axpy(1.0, t1); CHKERRQ(info);
  TaoFunctionReturn(0);
}

