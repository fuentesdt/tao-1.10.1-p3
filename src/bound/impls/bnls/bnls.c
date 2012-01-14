/*$Id: s.bnls.c 1.128 02/08/16 18:01:18-05:00 benson@rockies.mcs.anl.gov $*/

#include "bnls.h"       /*I "tao_solver.h" I*/



#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_BNLS"
static int TaoSolve_BNLS(TAO_SOLVER tao, void*solver){

  TAO_BNLS *bnls = (TAO_BNLS *)solver;
  int info;
  TaoInt lsflag,iter=0;
  TaoTerminateReason reason=TAO_CONTINUE_ITERATING;
  double f,f_full,gnorm,gdx,stepsize=1.0;
  TaoTruth success;
  TaoVec *XU, *XL;
  TaoVec *X,  *G=bnls->G, *PG=bnls->PG;
  TaoVec *R=bnls->R, *DXFree=bnls->DXFree;
  TaoVec *DX=bnls->DX, *Work=bnls->Work;
  TaoMat *H, *Hsub=bnls->Hsub;
  TaoIndexSet *FreeVariables = bnls->FreeVariables;

  TaoFunctionBegin;

  /* Check if upper bound greater than lower bound. */
  info = TaoGetSolution(tao,&X);CHKERRQ(info); bnls->X=X;
  info = TaoGetVariableBounds(tao,&XL,&XU);CHKERRQ(info);
  info = TaoEvaluateVariableBounds(tao,XL,XU); CHKERRQ(info);
  info = TaoGetHessian(tao,&H);CHKERRQ(info); bnls->H=H;

  /*   Project the current point onto the feasible set */
  info = X->Median(XL,X,XU); CHKERRQ(info);
  
  info = TaoComputeMeritFunctionGradient(tao,X,&f,G);CHKERRQ(info);
  
  while (reason==TAO_CONTINUE_ITERATING){
    
    /* Project the gradient and calculate the norm */
    info = PG->BoundGradientProjection(G,XL,X,XU);CHKERRQ(info);
    info = PG->Norm2(&gnorm); CHKERRQ(info);
    
    info = TaoMonitor(tao,iter++,f,gnorm,0.0,stepsize,&reason);
    CHKERRQ(info);
    if (reason!=TAO_CONTINUE_ITERATING) break;

    info = FreeVariables->WhichEqual(PG,G); CHKERRQ(info);

    info = TaoComputeHessian(tao,X,H);CHKERRQ(info);
    
    /* Create a reduced linear system */

    info = R->SetReducedVec(G,FreeVariables);CHKERRQ(info);
    info = R->Negate();CHKERRQ(info);

    info = DXFree->SetReducedVec(DX,FreeVariables);CHKERRQ(info);
    info = DXFree->SetToZero(); CHKERRQ(info);
    
    info = Hsub->SetReducedMatrix(H,FreeVariables,FreeVariables);CHKERRQ(info);

    bnls->gamma_factor /= 2;
    success = TAO_FALSE;

    while (success==TAO_FALSE) {
      
      /* Approximately solve the reduced linear system */
      info = TaoPreLinearSolve(tao,Hsub);CHKERRQ(info);
      info = TaoLinearSolve(tao,Hsub,R,DXFree,&success);CHKERRQ(info);

      info = DX->SetToZero(); CHKERRQ(info);
      info = DX->ReducedXPY(DXFree,FreeVariables);CHKERRQ(info);
      info = DX->Dot(G,&gdx); CHKERRQ(info);

      if (gdx>=0 || success==TAO_FALSE) { /* Modify diagonal of Hessian if not a descent direction */
        bnls->gamma_factor *= 2; 
        bnls->gamma = bnls->gamma_factor*(gnorm); 
#if !defined(PETSC_USE_COMPLEX)
        info=PetscInfo2(tao,"TaoSolve_NLS:  modify diagonal (assume same nonzero structure), gamma_factor=%g, gamma=%g\n",bnls->gamma_factor,bnls->gamma);
	CHKERRQ(info);
#else
        info=PetscInfo3(tao,"TaoSolve_NLS:  modify diagonal (asuume same nonzero structure), gamma_factor=%g, gamma=%g\n",
	     bnls->gamma_factor,PetscReal(bnls->gamma));CHKERRQ(info);
#endif
        info = Hsub->ShiftDiagonal(bnls->gamma);CHKERRQ(info);
	success = TAO_FALSE;
	
      } else {
	success = TAO_TRUE;
      }

    }
    
    stepsize=1.0;	
    info = TaoLineSearchApply(tao,X,G,DX,Work,
			      &f,&f_full,&stepsize,&lsflag);
    CHKERRQ(info);

    
  }  /* END MAIN LOOP  */

  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_BNLS"
static int TaoSetDown_BNLS(TAO_SOLVER tao, void*solver)
{
  TAO_BNLS *bnls = (TAO_BNLS *)solver;
  int      info;
  /* Free allocated memory in BNLS structure */
  TaoFunctionBegin;
  
  info = TaoVecDestroy(bnls->DX);CHKERRQ(info);bnls->DX=0;
  info = TaoVecDestroy(bnls->Work);CHKERRQ(info);
  info = TaoVecDestroy(bnls->DXFree);CHKERRQ(info);
  info = TaoVecDestroy(bnls->R);CHKERRQ(info);
  info = TaoVecDestroy(bnls->G);CHKERRQ(info);
  info = TaoVecDestroy(bnls->PG);CHKERRQ(info);
  info = TaoVecDestroy(bnls->XL);CHKERRQ(info);
  info = TaoVecDestroy(bnls->XU);CHKERRQ(info);
  
  info = TaoIndexSetDestroy(bnls->FreeVariables);CHKERRQ(info);
  info = TaoMatDestroy(bnls->Hsub);CHKERRQ(info);
  info = TaoDestroyLinearSolver(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_BNLS"
static int TaoSetOptions_BNLS(TAO_SOLVER tao, void*solver)
{
  int        info;
  TaoInt     ival;
  TaoTruth flg;

  TaoFunctionBegin;

  info = TaoOptionsHead("Newton Line Search Method for bound constrained optimization");CHKERRQ(info);

  info = TaoOptionInt("-redistribute","Redistribute Free variables (> 1 processors, only)","TaoPetscISType",1,&ival,&flg); CHKERRQ(info);

  info = TaoOptionName("-submatrixfree","Mask full matrix instead of extract submatrices","TaoPetscISType",&flg); CHKERRQ(info);

  info = TaoOptionsTail();CHKERRQ(info);
  info = TaoLineSearchSetFromOptions(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_BNLS"
static int TaoView_BNLS(TAO_SOLVER tao,void*solver)
{
  int        info;

  TaoFunctionBegin;
  info = TaoLineSearchView(tao);CHKERRQ(info);
  TaoFunctionReturn(0);
}


/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_BNLS"
static int TaoSetUp_BNLS(TAO_SOLVER tao, void*solver){

  int info;
  TAO_BNLS *bnls = (TAO_BNLS *)solver;
  TaoVec* X;
  TaoMat *HH;

  TaoFunctionBegin;
  info = TaoGetSolution(tao,&bnls->X);CHKERRQ(info); X=bnls->X;
  info = TaoGetHessian(tao,&bnls->H);CHKERRQ(info);  HH=bnls->H;

  /* Allocate some arrays */
  info = X->Clone(&bnls->DX); CHKERRQ(info);
  info = X->Clone(&bnls->Work); CHKERRQ(info);
  info = X->Clone(&bnls->DXFree); CHKERRQ(info);
  info = X->Clone(&bnls->R); CHKERRQ(info);
  info = X->Clone(&bnls->G); CHKERRQ(info);
  info = X->Clone(&bnls->PG); CHKERRQ(info);
  info = X->Clone(&bnls->XL); CHKERRQ(info);
  info = X->Clone(&bnls->XU); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao,bnls->PG);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,bnls->DX);CHKERRQ(info);
  info = TaoSetVariableBounds(tao,bnls->XL,bnls->XU);CHKERRQ(info);

  info = X->CreateIndexSet(&bnls->FreeVariables); CHKERRQ(info);
  info = bnls->H->CreateReducedMatrix(bnls->FreeVariables,bnls->FreeVariables,&bnls->Hsub); CHKERRQ(info);

  info = TaoCreateLinearSolver(tao,HH,100,0); CHKERRQ(info);

  info = TaoCheckFGH(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoGetDualVariables_BNLS" 
int TaoGetDualVariables_BNLS(TAO_SOLVER tao, TaoVec* DXL, TaoVec* DXU, void* solver)
{

  TAO_BNLS *bnls = (TAO_BNLS *) solver;
  TaoVec  *G=bnls->G,*GP=bnls->Work;
  int       info;

  TaoFunctionBegin;

  info = DXL->Waxpby(-1,G,1.0,GP); CHKERRQ(info);
  info = DXU->SetToZero(); CHKERRQ(info);
  info = DXL->PointwiseMaximum(DXL,DXU); CHKERRQ(info);

  info = DXU->Waxpby(-1.0,GP,1.0,G); CHKERRQ(info);
  info = DXU->Axpy(1.0,DXL); CHKERRQ(info);

  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_BNLS"
int TaoCreate_BNLS(TAO_SOLVER tao)
{
  TAO_BNLS *bnls;
  int      info;

  TaoFunctionBegin;

  info = TaoNew(TAO_BNLS,&bnls); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_BNLS)); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_BNLS,(void*)bnls); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_BNLS,TaoSetDown_BNLS); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetOptions_BNLS); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_BNLS); CHKERRQ(info);
  info=TaoSetTaoDualVariablesRoutine(tao,TaoGetDualVariables_BNLS); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,500); CHKERRQ(info);
  info = TaoSetTolerances(tao,1e-12,1e-12,0,0); CHKERRQ(info);

  /* Initialize pointers and variables */

  bnls->gamma = 0.0;
  bnls->gamma_factor = 0.01;
  bnls->DX=0;
  bnls->DXFree=0;
  bnls->R=0;
  bnls->Work=0;
  bnls->FreeVariables=0;
  bnls->Hsub=0;

  info = TaoCreateMoreThuenteBoundLineSearch(tao,0,0.9); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END
