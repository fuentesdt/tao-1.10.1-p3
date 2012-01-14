/*$Id$*/

#include "src/bound/impls/bcg/bcg.h"             /*I "tao_solver.h" I*/
char bcgtypename[]="tao_bcg";

/* ---------------------------------------------------------- */
#undef __FUNC__  
#define __FUNC__ "TaoSetUp_BCG"
int TaoSetUp_BCG(TAO_SOLVER tao,void *solver){

  TAO_BCG *cg = (TAO_BCG *) solver;
  TaoVec*    X;
  int    info;

  TaoFunctionBegin;
  
  info = TaoCheckFG(tao);CHKERRQ(info);
  info = TaoGetSolution(tao,&X);CHKERRQ(info);
    
  info = X->clone(&cg->Work); CHKERRQ(info);
  info = X->clone(&cg->DX); CHKERRQ(info);
  info = X->clone(&cg->Gprev); CHKERRQ(info);
  info = X->clone(&cg->GP); CHKERRQ(info);

  info = TaoLineSearchSetUp(tao);CHKERRQ(info);

  tao->vec_sol_update=cg->DX;

  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNC__  
#define __FUNC__ "TaoDestroy_BCG"
static int TaoDestroy_BCG(TAO_SOLVER tao, void *solver)
{
  TAO_BCG *cg = (TAO_BCG *) solver;
  int    info;

  TaoFunctionBegin;

  if (tao->setupcalled) {
    delete cg->Work;
    delete cg->DX;
    delete cg->Gprev;
    delete cg->GP;
  }
  info = TaoLineSearchDestroy(tao);CHKERRQ(info);
  TaoFree(cg);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/

#undef __FUNC__  
#define __FUNC__ "TaoSolve_BCG"
static int TaoSolve_BCG(TAO_SOLVER tao, void *cgptr)
{
  TAO_BCG *cg = (TAO_BCG *) cgptr;
  TaoVec*    X,*G; /* solution vector,  gradient vector */
  TaoVec*    Gprev=cg->Gprev, *GP=cg->GP;
  TaoVec*    DX=cg->DX, *Work=cg->Work;
  TaoVec  *XL,*XU;
  int    iter=0,lsflag=0,info;
  double gnorm2Prev,gdotgprev,gdx;
  double zero=0.0, minus_one = -1.0;
  double f_old,f,gnorm2,step=0;
  TaoTerminateReason reason;

  TaoFunctionBegin;
  info=TaoGetSolution(tao,&X);CHKERRQ(info);
  info=TaoGetGradient(tao,&G);CHKERRQ(info);

  info = TaoGetVariableBounds(tao,&XL,&XU); CHKERRQ(info);
  info = X->median(XL,X,XU); CHKERRQ(info);

  info = TaoComputeFunctionGradient(tao,X,&f,G);CHKERRQ(info);
  info = GP->boundGradientProjection(G,XL,X,XU); CHKERRQ(info);
  info = GP->norm2squared(&gnorm2);  CHKERRQ(info);

  info = DX->setToZero(); CHKERRQ(info); 
  info = Gprev->copyFrom(GP); CHKERRQ(info);

  cg->restarts=0;
  gnorm2Prev = gnorm2;

  /* Enter loop */
  while (1){

    /* Test for convergence */
    info = TaoMonitor(tao,iter++,f,gnorm2,0.0,step,&reason);CHKERRQ(info);
    if (reason!=TAO_CONTINUE_ITERATING) break;

    /* Determine beta, depending on method */
    info = GP->dot(Gprev,&gdotgprev); CHKERRQ(info);
    if (cg->type==TAO_CG_FletcherReeves){
      cg->beta=(gnorm2)/(gnorm2Prev);
    } else if (cg->type==TAO_CG_PolakRibiere){
      cg->beta=( (gnorm2)-gdotgprev )/(gnorm2Prev);
    } else {
      cg->beta=( (gnorm2)-gdotgprev )/(gnorm2Prev);
      if (cg->beta<0.0){
	cg->beta=0.0;
      }
    }

    /* Employ occasional restarts when successive gradients not orthogonal */
    if ( fabs(gdotgprev)/(gnorm2) > cg->eta || iter==0){ 
      printf("RESTART Beta: %4.2e\n",cg->beta);
      cg->beta=0.0;
    }

    if (cg->beta==0){
      cg->restarts++;
      PLogInfo(tao,"TaoCG: Restart CG at iterate %d with gradient direction.\n",tao->iter);
    }
    info = DX->scale(cg->beta); CHKERRQ(info);

    info = DX->negate(); CHKERRQ(info);
    info = DX->boundGradientProjection(DX,XL,X,XU); CHKERRQ(info);
    info = DX->negate(); CHKERRQ(info);

    info = DX->Axpy(minus_one,G); CHKERRQ(info);

    info = Gprev->copyFrom(GP); CHKERRQ(info);
    gnorm2Prev = gnorm2;

    info = Work->copyFrom(DX); CHKERRQ(info);
    info = Work->negate(); CHKERRQ(info);
    info = Work->boundGradientProjection(Work,XL,X,XU); CHKERRQ(info);
    info = Work->negate(); CHKERRQ(info);
    info = Work->dot(G,&gdx); CHKERRQ(info);
    if (cg->beta!=0 && gdx>=0){

      info = DX->copyFrom(GP); CHKERRQ(info);
      info = DX->negate(); CHKERRQ(info);
      cg->restarts++;
    } else {

    }
    info = DX->dot(G,&gdx); CHKERRQ(info);

    
    /* Line Search */
    step=1.5*step;
    step=TaoMax(1.5*step,0.1);
    step=1.0;

    info = TaoLineSearchApply(tao,X,G,DX,Work,&f,&step,&gdx,&lsflag);

    info = GP->boundGradientProjection(G,XL,X,XU); CHKERRQ(info);
    info = GP->norm2squared(&gnorm2); CHKERRQ(info);
    
  }
  
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNC__  
#define __FUNC__ "TaoSetOptions_BCG"
static int TaoSetOptions_BCG(TAO_SOLVER tao, void*solver)
{
  TAO_BCG *cg = (TAO_BCG *)solver;
  double   dtmp;
  int      info;
  char     method[20];
  PetscTruth flg;

  TaoFunctionBegin;

  info = OptionsHead("Nonlinear Conjugate Gradient method for unconstrained optimization");CHKERRQ(info);

  info = OptionsDouble("-tao_cg_restart","Restart method based on orthogonality of successive gradients","TaoBCGSetRestartTol",cg->eta,&dtmp, &flg);CHKERRQ(info);
  if (flg){
    info =TaoCGSetRestartTol(tao,dtmp);
  }

  info = OptionsString("-tao_cg_type [fletcher-reeves,polak-ribiere,polak-ribiere-plus]","Set Nonlinear BCG method","TaoBCGSetMethod","polak-ribiere-plus",method,19,&flg);
  CHKERRQ(info);
  if (flg){
    /*
    info = TaoBCGSetMethod(tao, method);CHKERRQ(info);
    */
  }
  info = TaoLineSearchSetFromOptions(tao);CHKERRQ(info);
  info = OptionsTail();CHKERRQ(info);

  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
#undef __FUNC__  
#define __FUNC__ "TaoView_BCG"
static int TaoView_BCG(TAO_SOLVER tao,void* solver)
{
  TAO_BCG   *cg = (TAO_BCG *)solver;
  int      info;

  TaoFunctionBegin;

  if (cg->type == TAO_CG_PolakRibiere){
    info = TaoPrintf(tao,"  tao_cg_type=polak-ribiere\n");CHKERRQ(info);      
  } else if (cg->type == TAO_CG_PRplus){
    info = TaoPrintf(tao,"  tao_cg_type=polak-ribiere-plus\n");CHKERRQ(info);
  } else if (cg->type == TAO_CG_FletcherReeves){
    info = TaoPrintf(tao,"  tao_cg_type=fletcher-reeves\n");CHKERRQ(info);
  }
  
  info = TaoPrintf2(tao,"  cg restarts=%d,   tao_cg_restart=%g\n",
		    cg->restarts,cg->eta);CHKERRQ(info);
  
  info = TaoLineSearchView(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
/*
    TAOCreate_BCG - Creates the data structure for the nonlinear BCG method
    and sets the function pointers for all the routines it needs to call
    (TAOSolve_BCG() etc.)

    It must be wrapped in EXTERN_C_BEGIN to be dynamically linkable in C++
*/
EXTERN_C_BEGIN
#undef __FUNC__  
#define __FUNC__ "TaoCreate_BCG"
int TaoCreate_BCG(TAO_SOLVER tao)
{
  TAO_BCG *cg;
  int    info;

  TaoFunctionBegin;

  cg = TaoNew(TAO_BCG);CHKPTRQ(cg);
  PLogObjectMemory(tao,sizeof(TAO_BCG));

  info = TaoSetSolver(tao,TaoSetUp_BCG,TaoSetOptions_BCG,TaoSolve_BCG,
		      TaoView_BCG,TaoDestroy_BCG,(void*)cg); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,2000); CHKERRQ(info);
  info = TaoSetTolerances(tao,1e-4,0,0,0); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,4000); CHKERRQ(info);

  info = TaoCreateProjectedLineSearch(tao); CHKERRQ(info);

  cg->eta = 100.0;
  cg->type = TAO_CG_PRplus;

  info = PetscObjectComposeFunctionDynamic((PetscObject)tao,"TaoBCGSetRestartTol_C",
					   "TaoBCGSetRestartTol_TaoBCG",
					   (void*)TaoBCGSetRestartTol_TaoBCG);CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END

EXTERN_C_BEGIN
#undef __FUNC__  
#define __FUNC__ "TaoCGSetRestartTol_TaoBCG" 
int TaoBCGSetRestartTol_TaoBCG(TAO_SOLVER tao,double eta)
{
  TAO_BCG *cg = (TAO_BCG *) tao->data;

  TaoFunctionBegin;
  cg->eta=TaoMax(eta,0.0);
  TaoFunctionReturn(0);
}
EXTERN_C_END

/* ---------------------------------------------------------- */
#undef __FUNC__  
#define __FUNC__ "TaoBCGSetRestartTol"
/*@ 
   TaoBCGSetRestartTol - Set the nonlinear conjugate gradient restart
   tolerance.  The algorithm restarts when the gradient at the current 
   point, g(x^k), and the gradient of the previous point, g^{k-1}, 
   satisfy the following inequality: 
   | g(x^k)^T g(x^{k-1}) | / g(x^k)^T g(x^k) > eta.

   Collective on TAO_SOLVER

   Input Parameters:
+  tao - the TAO_SOLVER context
-  eta - the number of previous points and gradients to be used.

   Output Parameter:

   Level: intermediate

   Options Database Keys: 
.  -tao_cg_restart <eta>

@*/
int TaoBCGSetRestartTol(TAO_SOLVER tao,double eta)
{
  int info,(*f)(TAO_SOLVER,double);

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE);
  if (eta < 0) SETERRQ(1,1,"Restart tolerance must be nonnegative");
  info = PetscObjectQueryFunction((PetscObject)tao,"TaoBCGSetRestartTol_C",(void **)&f);CHKERRQ(info);
  if (f){
    info = (*f)(tao,eta);CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "TaoBCGSetMethod" 
/* ---------------------------------------------------------- */
/*@
   TaoBCGSetMethod - Set the type of nonlinear conjugate gradient 
   method.  The three diffent types are "polak-ribiere",
   "polak-ribiere-plus", and "fletcher-reeves".
 
  Collective on TAO_SOLVER

   Input Parameters:
.  tao - the TAO_SOLVER context
.  method - either "polak-ribiere","polak-ribiere-plus", or "fletcher-reeves".

   Output Parameter:

   Options Database Keys: 
.  -tao_cg_type [ fletcher-reeves,polak-ribiere,polak-ribiere-plus ]


   Level: intermediate

@*/
int TaoBCGSetMethod(TAO_SOLVER tao,TAO_CGTYPES method)
{
  TAO_BCG *cg = (TAO_BCG *) tao->data;
  int info;
  TaoTruth issame;
  PetscTruth flg;

  TaoFunctionBegin;

  info = TaoCompareMethod(tao,bcgtypename,&issame);CHKERRQ(info);
  if (issame == TAO_FALSE) {
    TaoFunctionReturn(0);
  } else if (method==TAO_CG_PolakRibiere || method==TAO_CG_PRplus ||
	     method==TAO_CG_FletcherReeves) {
    cg->type=method;
  } else {
    SETERRQ(PETSC_ERR_ARG_WRONGSTATE,0,
	    "Error: TaoCGSetMethod: Unrecognized method %s\n");
  }

  TaoFunctionReturn(0);
}








