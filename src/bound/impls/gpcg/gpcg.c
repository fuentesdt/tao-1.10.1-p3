/*$Id$*/

#include "gpcg.h"        /*I "tao_solver.h" I*/

char gpcgname[]="tao_gpcg";
TaoMethod gpcgtypename = gpcgname;

static int TaoGradProjections(TAO_SOLVER, TAO_GPCG *);
static int GPCGCheckOptimalFace(TaoVec*, TaoVec*, TaoVec*, TaoVec*, TaoVec*, 
				TaoIndexSet*, TaoIndexSet*, TaoTruth *optimal);

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_GPCG"
static int TaoSetDown_GPCG(TAO_SOLVER tao, void*solver)
{
  TAO_GPCG *gpcg = (TAO_GPCG *)solver;
  int      info;
  /* Free allocated memory in GPCG structure */
  TaoFunctionBegin;
  
  info = TaoVecDestroy(gpcg->X_New);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->G_New);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->DX);CHKERRQ(info);gpcg->DX=0;
  info = TaoVecDestroy(gpcg->Work);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->DXFree);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->R);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->B);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->G);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->PG);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->XL);CHKERRQ(info);
  info = TaoVecDestroy(gpcg->XU);CHKERRQ(info);
  
  info = TaoIndexSetDestroy(gpcg->Free_Local);CHKERRQ(info);
  info = TaoIndexSetDestroy(gpcg->TT);CHKERRQ(info);
  info = TaoMatDestroy(gpcg->Hsub);CHKERRQ(info);
  
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetFromOptions_GPCG"
static int TaoSetFromOptions_GPCG(TAO_SOLVER tao, void*solver)
{
  TAO_GPCG *gpcg = (TAO_GPCG *)solver;
  int      info;
  TaoInt   ival;
  TaoTruth flg;

  TaoFunctionBegin;
  info = TaoOptionsHead("Gradient Projection, Conjugate Gradient method for bound constrained optimization");CHKERRQ(info);

  info=TaoOptionInt("-gpcg_maxpgits","maximum number of gradient projections per GPCG iterate",0,gpcg->maxgpits,&gpcg->maxgpits,&flg);
  CHKERRQ(info);

  info = TaoOptionInt("-redistribute","Redistribute Free variables (> 1 processors, only)","TaoPetscISType",1,&ival,&flg); CHKERRQ(info);

  info = TaoOptionName("-submatrixfree","Mask full matrix instead of extract submatrices","TaoPetscISType",&flg); CHKERRQ(info);


  info = TaoOptionsTail();CHKERRQ(info);
  info=TaoLineSearchSetFromOptions(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_GPCG"
static int TaoView_GPCG(TAO_SOLVER tao,void*solver)
{
  TAO_GPCG *gpcg = (TAO_GPCG *)solver;
  int      info;

  TaoFunctionBegin;

  info = TaoPrintInt(tao," Total PG its: %d,",gpcg->total_gp_its);CHKERRQ(info);
  info = TaoPrintDouble(tao," PG tolerance: %4.3f \n",gpcg->pg_ftol);CHKERRQ(info);
  info = TaoLineSearchView(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoGPCGComputeFunctionGradient"
int TaoGPCGComputeFunctionGradient(TAO_SOLVER tao, TaoVec *XX, double *ff, TaoVec *GG){
  TAO_GPCG *gpcg;
  TaoMat *HH;
  int info;
  double f1,f2;

  TaoFunctionBegin;
  info = TaoGetSolverContext(tao,gpcgtypename,(void**)&gpcg); CHKERRQ(info);
  if (gpcg==0){TaoFunctionReturn(0);}
  info = TaoGetHessian(tao,&HH);CHKERRQ(info);
  info = HH->Multiply(XX,GG);CHKERRQ(info);
  info = XX->Dot(GG,&f1);CHKERRQ(info);
  info = XX->Dot(gpcg->B,&f2);CHKERRQ(info);
  info = GG->Axpy(1.0,gpcg->B);CHKERRQ(info);
  *ff=f1/2.0 + f2 + gpcg->c;
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_GPCG"
static int TaoSetUp_GPCG(TAO_SOLVER tao,void*solver){

  int      info;
  TaoInt   n;
  TAO_GPCG *gpcg = (TAO_GPCG *) solver;
  TaoVec   *X;
  TaoMat   *HH;
  TaoIndexSet *TIS;

  TaoFunctionBegin;

  info = TaoGetSolution(tao,&X);CHKERRQ(info); gpcg->X=X;
  info = TaoGetHessian(tao,&HH);CHKERRQ(info); gpcg->H=HH;

  /* Allocate some arrays */
  info=X->Clone(&gpcg->DX); CHKERRQ(info);
  info=X->Clone(&gpcg->B); CHKERRQ(info);
  info=X->Clone(&gpcg->Work); CHKERRQ(info);
  info=X->Clone(&gpcg->X_New); CHKERRQ(info);
  info=X->Clone(&gpcg->G_New); CHKERRQ(info);
  info=X->Clone(&gpcg->DXFree); CHKERRQ(info);
  info=X->Clone(&gpcg->R); CHKERRQ(info);
  info=X->Clone(&gpcg->G); CHKERRQ(info);
  info=X->Clone(&gpcg->PG); CHKERRQ(info);
  info=X->Clone(&gpcg->XL); CHKERRQ(info);
  info=X->Clone(&gpcg->XU); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao,gpcg->PG);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,gpcg->DX);CHKERRQ(info);
  info = TaoSetVariableBounds(tao,gpcg->XL,gpcg->XU);CHKERRQ(info);

  info = X->GetDimension(&n); CHKERRQ(info);
  gpcg->n=n;
  info = TaoCreateLinearSolver(tao,HH,300,0); CHKERRQ(info);

  info = X->CreateIndexSet(&TIS); CHKERRQ(info);
  gpcg->Free_Local = TIS;
  info = gpcg->Free_Local->Duplicate(&gpcg->TT); CHKERRQ(info);

  info = HH->CreateReducedMatrix(TIS,TIS,&gpcg->Hsub); CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_GPCG"
static int TaoSolve_GPCG(TAO_SOLVER tao, void *solver)
{
  TAO_GPCG *gpcg = (TAO_GPCG *)solver ;
  int info;
  TaoInt lsflag,iter=0;
  TaoTruth optimal_face=TAO_FALSE,success;
  double actred,f,f_new,f_full,gnorm,gdx,stepsize;
  double c;
  TaoVec *XU, *XL;
  TaoVec *X,  *G=gpcg->G , *B=gpcg->B, *PG=gpcg->PG;
  TaoVec *R=gpcg->R, *DXFree=gpcg->DXFree;
  TaoVec *G_New=gpcg->G_New;
  TaoVec *DX=gpcg->DX, *Work=gpcg->Work;
  TaoMat *H, *Hsub=gpcg->Hsub;
  TaoIndexSet *Free_Local = gpcg->Free_Local, *TIS=gpcg->TT;
  TaoTerminateReason reason;

  TaoFunctionBegin;

  /* Check if upper bound greater than lower bound. */
  info = TaoGetSolution(tao,&X);CHKERRQ(info);
  info = TaoGetHessian(tao,&H);CHKERRQ(info);

  info = TaoGetVariableBounds(tao,&XL,&XU);CHKERRQ(info);
  info = TaoEvaluateVariableBounds(tao,XL,XU); CHKERRQ(info);
  info = X->Median(XL,X,XU); CHKERRQ(info);

  info = TaoComputeHessian(tao,X,H); CHKERRQ(info);
  info = TaoComputeFunctionGradient(tao,X,&f,B);
  CHKERRQ(info);

  /* Compute quadratic representation */
  info = H->Multiply(X,Work); CHKERRQ(info);
  info = X->Dot(Work,&c); CHKERRQ(info);
  info = B->Axpy(-1.0,Work); CHKERRQ(info);
  info = X->Dot(B,&stepsize); CHKERRQ(info);
  gpcg->c=f-c/2.0-stepsize;

  info = Free_Local->WhichBetween(XL,X,XU); CHKERRQ(info);
  
  info = TaoGPCGComputeFunctionGradient(tao, X, &gpcg->f , G);
  
  /* Project the gradient and calculate the norm */
  info = G_New->CopyFrom(G);CHKERRQ(info);
  info = PG->BoundGradientProjection(G,XL,X,XU);CHKERRQ(info);
  info = PG->Norm2(&gpcg->gnorm); CHKERRQ(info);
  gpcg->step=1.0;

    /* Check Stopping Condition      */
  info=TaoMonitor(tao,iter++,gpcg->f,gpcg->gnorm,0,gpcg->step,&reason); CHKERRQ(info);

  while (reason == TAO_CONTINUE_ITERATING){

    info = TaoGradProjections(tao, gpcg); CHKERRQ(info);

    info = Free_Local->WhichBetween(XL,X,XU); CHKERRQ(info);
    info = Free_Local->GetSize(&gpcg->n_free); CHKERRQ(info);
    f=gpcg->f; gnorm=gpcg->gnorm; 

    if (gpcg->n_free > 0){
      
      /* Create a reduced linear system */
      info = R->SetReducedVec(G,Free_Local);CHKERRQ(info);
      info = R->Negate(); CHKERRQ(info);
      info = DXFree->SetReducedVec(DX,Free_Local);CHKERRQ(info);
      info = DXFree->SetToZero(); CHKERRQ(info);

      info = Hsub->SetReducedMatrix(H,Free_Local,Free_Local);CHKERRQ(info);

      info = TaoPreLinearSolve(tao,Hsub);CHKERRQ(info);

      /* Approximately solve the reduced linear system */
      info = TaoLinearSolve(tao,Hsub,R,DXFree,&success);CHKERRQ(info);
      
      info=DX->SetToZero(); CHKERRQ(info);
      info=DX->ReducedXPY(DXFree,Free_Local);CHKERRQ(info);
      
      info = G->Dot(DX,&gdx); CHKERRQ(info);
      
      stepsize=1.0; f_new=f;
      info = TaoLineSearchApply(tao,X,G,DX,Work,
				&f_new,&f_full,&stepsize,&lsflag);
      CHKERRQ(info);
      
      actred = f_new - f;
      
      /* Evaluate the function and gradient at the new point */      
      info =  PG->BoundGradientProjection(G,XL,X,XU);
      CHKERRQ(info);
      info = PG->Norm2(&gnorm);  CHKERRQ(info);      
      f=f_new;
      
      info = GPCGCheckOptimalFace(X,XL,XU,PG,Work, Free_Local, TIS,
				  &optimal_face); CHKERRQ(info);
      
    } else {
      
      actred = 0; stepsize=1.0;
      /* if there were no free variables, no cg method */

    }

    info = TaoMonitor(tao,iter,f,gnorm,0.0,stepsize,&reason); CHKERRQ(info);
    gpcg->f=f;gpcg->gnorm=gnorm; gpcg->actred=actred;
    if (reason!=TAO_CONTINUE_ITERATING) break;
    iter++;


  }  /* END MAIN LOOP  */

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGradProjections"
static int TaoGradProjections(TAO_SOLVER tao, TAO_GPCG *gpcg)
{
  int info;
  TaoInt lsflag=0,i;
  TaoTruth optimal_face=TAO_FALSE;
  double actred=-1.0,actred_max=0.0, gAg,gtg=gpcg->gnorm,alpha;
  double f_new,f_full,gdx;
  TaoMat *H;
  TaoVec *DX=gpcg->DX,*XL=gpcg->XL,*XU=gpcg->XU,*Work=gpcg->Work;
  TaoVec *X=gpcg->X,*G=gpcg->G;
  /*
     The gradient and function value passed into and out of this
     routine should be current and correct.
     
     The free, active, and binding variables should be already identified
  */
  
  TaoFunctionBegin;
  
  info = TaoGetSolution(tao,&X);CHKERRQ(info);
  info = TaoGetHessian(tao,&H);CHKERRQ(info);
  info = TaoGetVariableBounds(tao,&XL,&XU);CHKERRQ(info);

  for (i=0;i<gpcg->maxgpits;i++){

    if ( -actred <= (gpcg->pg_ftol)*actred_max) break;
 
    info = DX->BoundGradientProjection(G,XL,X,XU); CHKERRQ(info);
    info = DX->Negate(); CHKERRQ(info);
    info = DX->Dot(G,&gdx); CHKERRQ(info);

    info= H->Multiply(DX,Work); CHKERRQ(info);
    info= DX->Dot(Work,&gAg); CHKERRQ(info);
 
    gpcg->gp_iterates++; gpcg->total_gp_its++;    
  
    gtg=-gdx;
    alpha = TaoAbsDouble(gtg/gAg);
    gpcg->stepsize = alpha; f_new=gpcg->f;

    info = TaoLineSearchApply(tao,X,G,DX,Work,
			      &f_new,&f_full,&gpcg->stepsize,&lsflag);
    CHKERRQ(info);

    /* Update the iterate */
    actred = f_new - gpcg->f;
    actred_max = TaoMax(actred_max,-(f_new - gpcg->f));
    gpcg->f = f_new;
    info = GPCGCheckOptimalFace(X,XL,XU,G,Work,gpcg->Free_Local,gpcg->TT,
				&optimal_face); CHKERRQ(info);

    if ( optimal_face == TAO_TRUE ) break;

  }
  
  gpcg->gnorm=gtg;
  TaoFunctionReturn(0);

} /* End gradient projections */


#undef __FUNCT__  
#define __FUNCT__ "GPCGCheckOptimalFace"
static int GPCGCheckOptimalFace(TaoVec *X,TaoVec *XL,TaoVec*XU,TaoVec *PG,TaoVec*W,
				TaoIndexSet*Free_Local, TaoIndexSet*TT,
				TaoTruth *optimal)
{
  int info;
  TaoInt n_free;
  double rr;
  TaoTruth same;

  TaoFunctionBegin;
  *optimal = TAO_FALSE;

  /* Also check to see if the active set is the same */

  info = TT->WhichBetween(XL,X,XU); CHKERRQ(info);
  info = Free_Local->IsSame(TT,&same); CHKERRQ(info);
  info = Free_Local->GetSize(&n_free); CHKERRQ(info);
  if (same == TAO_FALSE){
    info = Free_Local->WhichBetween(XL,X,XU); CHKERRQ(info);
    *optimal = TAO_FALSE;
    TaoFunctionReturn(0);
  } else {
    *optimal = TAO_TRUE;
  }


  info = W->CopyFrom(PG); CHKERRQ(info);
  info = W->Negate(); CHKERRQ(info);

  info = W->BoundGradientProjection(W,XL,X,XU); CHKERRQ(info);
  info = W->Axpy(1.0,PG); CHKERRQ(info);

  info = W->Norm2(&rr); CHKERRQ(info);
  if (rr>0) *optimal = TAO_FALSE;

  *optimal = TAO_FALSE;
  /*
  info = gpcg->TT->whichNonNegative(W); CHKERRQ(info);
  info = gpcg->TT->GetSize(&n); CHKERRQ(info);
  if (n==0) *optimal = TAO_TRUE;
  */
  TaoFunctionReturn(0);
}

int TaoDefaultMonitor_GPCG(TAO_SOLVER tao,void *dummy)
{
  TAO_GPCG *gpcg;
  double   fct,gnorm;
  int      info;
  TaoInt   its,nfree;

  TaoFunctionBegin;
  info = TaoGetSolutionStatus(tao,&its,&fct,&gnorm,0,0,0);CHKERRQ(info);
  info = TaoGetSolverContext(tao,gpcgtypename,(void**)&gpcg); CHKERRQ(info);
  if (gpcg==0){TaoFunctionReturn(0);}
  nfree=gpcg->n_free;
  info = TaoPrintInt(tao,"iter: %d,",its);CHKERRQ(info);
  info = TaoPrintDouble(tao," Fcn value: %g,",fct);CHKERRQ(info);
  info = TaoPrintDouble(tao," PGrad. norm: %g, ",gnorm);CHKERRQ(info);
  info = TaoPrintInt(tao,"free vars:%d \n",nfree);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetDualVariables_GPCG" 
int TaoGetDualVariables_GPCG(TAO_SOLVER tao, TaoVec* DXL, TaoVec* DXU, void* solver)
{

  TAO_GPCG *gpcg = (TAO_GPCG *) solver;
  TaoVec  *G=gpcg->G,*GP=gpcg->Work;
  TaoVec  *X,*XL,*XU;
  int       info;

  TaoFunctionBegin;
  info = TaoGetVariableBounds(tao,&XL,&XU); CHKERRQ(info);
  info = TaoGetSolution(tao,&X); CHKERRQ(info);
  info = GP->BoundGradientProjection(G,XL,X,XU); CHKERRQ(info);

  info = DXL->Waxpby(-1.0,G,1.0,GP); CHKERRQ(info);
  info = DXU->SetToZero(); CHKERRQ(info);
  info = DXL->PointwiseMaximum(DXL,DXU); CHKERRQ(info);

  info = DXU->Waxpby(-1.0,GP,1.0,G); CHKERRQ(info);
  info = GP->SetToZero(); CHKERRQ(info);
  info = DXU->PointwiseMinimum(GP,DXU); CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_GPCG"
int TaoCreate_GPCG(TAO_SOLVER tao)
{
  TAO_GPCG *gpcg;
  int      info;

  TaoFunctionBegin;

  info = TaoNew(TAO_GPCG,&gpcg); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_GPCG)); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_GPCG,(void*)gpcg); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_GPCG,TaoSetDown_GPCG); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetFromOptions_GPCG); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_GPCG); CHKERRQ(info);
  info=TaoSetTaoDualVariablesRoutine(tao,TaoGetDualVariables_GPCG); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,500); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,100000); CHKERRQ(info);
  info = TaoSetTolerances(tao,1e-12,1e-12,0,0); CHKERRQ(info);

  /* Initialize pointers and variables */
  gpcg->n=0;
  gpcg->maxgpits = 8;
  gpcg->pg_ftol = 0.1;

  gpcg->gp_iterates=0; /* Cumulative number */
  gpcg->total_gp_its = 0;
 
  /* Initialize pointers and variables */
  gpcg->n_bind=0;
  gpcg->n_free = 0;
  gpcg->n_upper=0;
  gpcg->n_lower=0;

  //  info = TaoCreateProjectedLineSearch(tao); CHKERRQ(info);
  info = TaoGPCGCreateLineSearch(tao); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END




