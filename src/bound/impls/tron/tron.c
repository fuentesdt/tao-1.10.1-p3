/*$Id$*/

#include "tron.h"       /*I "tao_solver.h" I*/
#include "petscksp.h"
#include "petscpc.h"
#include "src/petsctao/linearsolver/taolinearsolver_petsc.h"
#include "src/petsctao/vector/taovec_petsc.h"
#include "private/kspimpl.h"
#include "private/pcimpl.h"

/* TRON Routines */
static int TaoGradProjections(TAO_SOLVER,TAO_TRON *);
static int TronCheckOptimalFace(TaoVec*, TaoVec*, TaoVec*, TaoVec*, TaoVec*, 
				TaoIndexSet*, TaoIndexSet*, TaoTruth *optimal);

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_TRON"
static int TaoSetDown_TRON(TAO_SOLVER tao, void*solver)
{
  TAO_TRON *tron = (TAO_TRON *)solver;
  int      info;
  /* Free allocated memory in TRON structure */
  TaoFunctionBegin;
  
  info = TaoVecDestroy(tron->X_New);CHKERRQ(info);
  info = TaoVecDestroy(tron->G_New);CHKERRQ(info);
  info = TaoVecDestroy(tron->DX);CHKERRQ(info);tron->DX=0;
  info = TaoVecDestroy(tron->Work);CHKERRQ(info);
  info = TaoVecDestroy(tron->DXFree);CHKERRQ(info);
  info = TaoVecDestroy(tron->R);CHKERRQ(info);
  info = TaoVecDestroy(tron->G);CHKERRQ(info);
  info = TaoVecDestroy(tron->PG);CHKERRQ(info);
  info = TaoVecDestroy(tron->XL);CHKERRQ(info);
  info = TaoVecDestroy(tron->XU);CHKERRQ(info);
  
  info = TaoIndexSetDestroy(tron->Free_Local);CHKERRQ(info);
  info = TaoIndexSetDestroy(tron->TT);CHKERRQ(info);
  info = TaoMatDestroy(tron->Hsub);CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao,0);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,0);CHKERRQ(info);
  info = TaoSetVariableBounds(tao,0,0);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_TRON"
static int TaoSetOptions_TRON(TAO_SOLVER tao, void*solver)
{
  TAO_TRON  *tron = (TAO_TRON *)solver;
  int        info;
  TaoInt     ival;
  TaoTruth flg;

  TaoFunctionBegin;

  info = TaoOptionsHead("Newton Trust Region Method for bound constrained optimization");CHKERRQ(info);

  info = TaoOptionInt("-tron_maxgpits","maximum number of gradient projections per TRON iterate","TaoSetMaxGPIts",tron->maxgpits,&tron->maxgpits,&flg);
  CHKERRQ(info);

  info = TaoOptionInt("-redistribute","Redistribute Free variables (> 1 processors, only)","TaoPetscISType",1,&ival,&flg); CHKERRQ(info);

  info = TaoOptionName("-submatrixfree","Mask full matrix instead of extract submatrices","TaoPetscISType",&flg); CHKERRQ(info);

  info = TaoOptionsTail();CHKERRQ(info);
  info = TaoLineSearchSetFromOptions(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_TRON"
static int TaoView_TRON(TAO_SOLVER tao,void*solver)
{
  TAO_TRON  *tron = (TAO_TRON *)solver;
  int        info;

  TaoFunctionBegin;
  /*
  info = TaoPrintf1(tao," Variables, Total: %d,",tron->n);
  info = TaoPrintf3(tao,"Free: %d,  Binding: %d \n",
		    tron->n_free, tron->n - tron->n_free,
		    tron->n_bind);CHKERRQ(info);
  info = TaoPrintf1(tao,"            Equal lower bound: %d,",
		    tron->n_lower);CHKERRQ(info);
  info = TaoPrintf1(tao,"  Equal upper bound: %d \n",
		    tron->n_upper);CHKERRQ(info);
  */
  info = TaoPrintInt(tao," Total PG its: %d,",tron->total_gp_its);CHKERRQ(info);
  info = TaoPrintDouble(tao," PG tolerance: %4.3f \n",tron->pg_ftol);CHKERRQ(info);
  info = TaoLineSearchView(tao);CHKERRQ(info);
  info = TaoPrintStatement(tao,"  Linear Solver minimizes quadratic over Trust Region: \n");CHKERRQ(info);

  TaoFunctionReturn(0);
}


/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_TRON"
static int TaoSetUp_TRON(TAO_SOLVER tao, void*solver){

  int info;
  TAO_TRON *tron = (TAO_TRON *)solver;
  TaoVec* X;
  TaoMat *HH;
  TaoIndexSet *TIS;

  TaoFunctionBegin;
  info = TaoGetSolution(tao,&tron->X);CHKERRQ(info); X=tron->X;
  info = TaoGetHessian(tao,&tron->H);CHKERRQ(info);  HH=tron->H;

  /* Allocate some arrays */
  info = X->Clone(&tron->DX); CHKERRQ(info);
  info = X->Clone(&tron->X_New); CHKERRQ(info);
  info = X->Clone(&tron->G_New); CHKERRQ(info);
  info = X->Clone(&tron->Work); CHKERRQ(info);
  info = X->Clone(&tron->DXFree); CHKERRQ(info);
  info = X->Clone(&tron->R); CHKERRQ(info);
  info = X->Clone(&tron->G); CHKERRQ(info);
  info = X->Clone(&tron->PG); CHKERRQ(info);
  info = X->Clone(&tron->XL); CHKERRQ(info);
  info = X->Clone(&tron->XU); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao,tron->PG);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,tron->DX);CHKERRQ(info);
  info = TaoSetVariableBounds(tao,tron->XL,tron->XU);CHKERRQ(info);

  info = X->GetDimension(&tron->n); CHKERRQ(info);
  
  info = X->CreateIndexSet(&tron->Free_Local); CHKERRQ(info);
  info = tron->Free_Local->Duplicate(&tron->TT); CHKERRQ(info);

  TIS=tron->Free_Local;
  info = tron->H->CreateReducedMatrix(TIS,TIS,&tron->Hsub); CHKERRQ(info);
  info = TaoCreateLinearSolver(tao,HH,220,0); CHKERRQ(info);

  info = TaoCheckFGH(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_TRON"
static int TaoSolve_TRON(TAO_SOLVER tao, void*solver){

  TAO_TRON *tron = (TAO_TRON *)solver;
  int info;
  TaoInt lsflag,iter=0;
  TaoTerminateReason reason;
  TaoTruth optimal_face=TAO_FALSE,success;
  double prered,actred,delta,f,f_new,f_full,rhok,gnorm,gdx,xdiff,stepsize;
  TaoVec *XU, *XL;
  TaoVec *X,  *G;
  TaoVec *PG=tron->PG;
  TaoVec *R=tron->R, *DXFree=tron->DXFree;
  TaoVec *X_New=tron->X_New, *G_New=tron->G_New;
  TaoVec *DX=tron->DX, *Work=tron->Work;
  TaoMat *H, *Hsub=tron->Hsub;
  TaoIndexSet *Free_Local = tron->Free_Local, *TIS=tron->TT;

  TaoFunctionBegin;

  // Get initial trust region radius
  info = TaoGetInitialTrustRegionRadius(tao, &tron->delta); CHKERRQ(info);
  if (tron->delta <= 0) {
    SETERRQ(1, "Initial trust region radius must be positive");
  }

  // Get vectors we will need
  info = TaoGetSolution(tao, &X); CHKERRQ(info);
  info = TaoGetGradient(tao, &G); CHKERRQ(info);
  info = TaoGetHessian(tao, &H); CHKERRQ(info);
  info = TaoGetVariableBounds(tao, &XL, &XU); CHKERRQ(info);

  // Check that upper bound greater than lower bound
  info = TaoEvaluateVariableBounds(tao, XL, XU); CHKERRQ(info);

  tron->pgstepsize=1.0;

  /*   Project the current point onto the feasible set */
  info = X->Median(XL,X,XU); CHKERRQ(info);
  
  info = TaoComputeMeritFunctionGradient(tao,X,&tron->f,G);CHKERRQ(info);
  info = Free_Local->WhichBetween(XL,X,XU); CHKERRQ(info);
  
  /* Project the gradient and calculate the norm */
  //  info = G_New->CopyFrom(G);CHKERRQ(info);
  info = PG->BoundGradientProjection(G,XL,X,XU);CHKERRQ(info);
  info = PG->Norm2(&tron->gnorm); CHKERRQ(info);

  if (tron->delta <= 0.0){
    tron->delta=TaoMax(tron->gnorm*tron->gnorm,1.0);
    //    tron->delta = TAO_INFINITY;
  }

  tron->stepsize=tron->delta;

  info=PetscInfo1(tao,"delta =%22.12e\n",tron->delta);

  info = TaoMonitor(tao,iter++,tron->f,tron->gnorm,0.0,tron->delta,&reason);
  CHKERRQ(info);

  while (reason==TAO_CONTINUE_ITERATING){
    
    info = TaoGradProjections(tao,tron); CHKERRQ(info);

    info = Free_Local->WhichBetween(XL,X,XU); CHKERRQ(info);
    info = Free_Local->GetSize(&tron->n_free); CHKERRQ(info);
    f=tron->f; delta=tron->delta; gnorm=tron->gnorm; 
    info=PetscInfo4(tao,"n_free=%d,f =%22.12e , delta =%22.12e , gnorm =%22.12e\n",tron->n_free,f,delta,gnorm);
    

    if (tron->n_free > 0){
      
      /* view and Modify the linear solver */
      TaoLinearSolver *tls;
      info = TaoGetLinearSolver(tao, &tls); CHKERRQ(info);
      TaoLinearSolverPetsc *pls;
      pls  = dynamic_cast <TaoLinearSolverPetsc *> (tls);
      KSP pksp = pls->GetKSP();
      PetscScalar ewAtol  = PetscMin(0.5,gnorm)*gnorm;
      info = KSPSetTolerances(pksp,PETSC_DEFAULT,ewAtol,
                      PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(info);
      pksp->printreason = PETSC_TRUE;
      info = KSPView(pksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(info);

      /* Compute the Hessian */
      info = TaoComputeHessian(tao,X,H);CHKERRQ(info);

      /* Create a reduced linear system */
      info = R->SetReducedVec(G,Free_Local);CHKERRQ(info);
      info = R->Negate(); CHKERRQ(info);
      info = DXFree->SetReducedVec(DX,Free_Local);CHKERRQ(info);
      info = DXFree->SetToZero(); CHKERRQ(info);

      info = Hsub->SetReducedMatrix(H,Free_Local,Free_Local);CHKERRQ(info);

      info = TaoPreLinearSolve(tao,Hsub);CHKERRQ(info);

      while (1) {

 	/* Approximately solve the reduced linear system */

	info = TaoLinearSolveTrustRegion(tao, Hsub, R, DXFree, delta, &success); CHKERRQ(info);

	info=DX->SetToZero(); CHKERRQ(info);
	info=DX->ReducedXPY(DXFree,Free_Local);CHKERRQ(info);

	info = G->Dot(DX,&gdx); CHKERRQ(info);
	info = PetscInfo1(tao,"Expected decrease in function value: %14.12e\n",gdx); CHKERRQ(info);

	stepsize=1.0; f_new=f;
	info = X_New->CopyFrom(X); CHKERRQ(info);
	info = G_New->CopyFrom(G); CHKERRQ(info);
	
	info = TaoLineSearchApply(tao,X_New,G_New,DX,Work,
				  &f_new,&f_full,&stepsize,&lsflag);
	CHKERRQ(info);
	info = H->Multiply(DX,Work); CHKERRQ(info);
	info = Work->Aypx(0.5,G); CHKERRQ(info);
	info = Work->Dot(DX,&prered); CHKERRQ(info);
	actred = f_new - f;
	
	if (actred<0) rhok=TaoAbsScalar(-actred/prered);
	else rhok=0.0;

	/* Compare actual improvement to the quadratic model */
	info = PetscInfo2(tao,"rhok=%14.12e delta=%14.12e\n",rhok,delta); CHKERRQ(info);
	if (rhok > tron->eta1) { /* Accept the point */

	  info = DX->Waxpby(1.0,X_New,-1.0, X); CHKERRQ(info);
	  info = DX->Norm2(&xdiff); CHKERRQ(info);
	  xdiff*=stepsize;

	  /* Adjust trust region size */
	  if (rhok < tron->eta2 ){
	    delta = TaoMin(xdiff,delta)*tron->sigma1;
	  } else if (rhok > tron->eta4 ){
	    delta= TaoMin(xdiff,delta)*tron->sigma3;
	  } else if (rhok > tron->eta3 ){
	    delta=TaoMin(xdiff,delta)*tron->sigma2;
	  }
	  info = PetscInfo1(tao,"update delta=%14.12e\n",delta); CHKERRQ(info);

	  info =  PG->BoundGradientProjection(G_New,XL,X_New,XU);
	  CHKERRQ(info);
	  info = PG->Norm2(&gnorm);  CHKERRQ(info);
	  info = TronCheckOptimalFace(X_New,XL,XU,G_New,PG, Free_Local, TIS,
				      &optimal_face); CHKERRQ(info);	  
          if (stepsize < 1 || optimal_face==TAO_FALSE || reason!=TAO_CONTINUE_ITERATING ){
            f=f_new;
            info = X->CopyFrom(X_New); CHKERRQ(info);
            info = G->CopyFrom(G_New); CHKERRQ(info);
            break;
          }
	  if (delta<=1e-30){
            break;
	  }
	} 
	else if (delta <= 1e-30) {
	  break;
	}
        else {
	  delta /= 4.0;
	}
      } /* end linear solve loop */
      
    } else {
      
      actred=0;
      info =  Work->BoundGradientProjection(G,XL,X,XU);
      CHKERRQ(info);
      info = Work->Norm2(&gnorm);  CHKERRQ(info);
      info = PetscInfo1(tao,"no free variables gnorm: %14.12e\n",gnorm); CHKERRQ(info);
      /* if there were no free variables, no cg method */

    }

    tron->f=f;tron->gnorm=gnorm; tron->actred=actred; tron->delta=delta;
    info = TaoMonitor(tao,iter,f,gnorm,0.0,delta,&reason); CHKERRQ(info);
    if (reason!=TAO_CONTINUE_ITERATING) break;
    iter++;
    
  }  /* END MAIN LOOP  */

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGradProjections"
static int TaoGradProjections(TAO_SOLVER tao,TAO_TRON *tron)
{
  int info;
  TaoInt lsflag=0,i;
  TaoTruth sameface=TAO_FALSE;
  double actred=-1.0,actred_max=0.0;
  double f_new, f_full;
  TaoVec *DX=tron->DX,*XL=tron->XL,*XU=tron->XU,*Work=tron->Work;
  TaoVec *X=tron->X,*G=tron->G;
  TaoIndexSet *TT1=tron->Free_Local, *TT2=tron->TT, *TT3;
  /*
     The gradient and function value passed into and out of this
     routine should be current and correct.
     
     The free, active, and binding variables should be already identified
  */
  
  TaoFunctionBegin;
  
  info = TaoGetSolution(tao,&X);CHKERRQ(info);
  info = TaoGetGradient(tao,&G);CHKERRQ(info);
  info = TaoGetVariableBounds(tao,&XL,&XU);CHKERRQ(info);

  info = TT1->WhichBetween(XL,X,XU); CHKERRQ(info);

  for (i=0;i<tron->maxgpits;i++){

    info = PetscInfo4(tao,"GradProjection %d, actred=%22.12e,  pg_ftol=%22.12e, actred_max=%22.12e \n",
                                           i,actred,tron->pg_ftol,actred_max); CHKERRQ(info);
    if ( -actred <= (tron->pg_ftol)*actred_max) break;
  
    tron->gp_iterates++; tron->total_gp_its++;      
    f_new=tron->f;

    info = DX->ScaleCopyFrom(-1.0,G); CHKERRQ(info);

    info = TaoLineSearchApply(tao,X,G,DX,Work,
			      &f_new,&f_full,&tron->pgstepsize,&lsflag);
    CHKERRQ(info);

    /* Update the iterate */
    actred = f_new - tron->f;
    actred_max = TaoMax(actred_max,-(f_new - tron->f));
    tron->f = f_new;

    info = TT2->WhichBetween(XL,X,XU); CHKERRQ(info);
    info = TT2->IsSame(TT1,&sameface);  CHKERRQ(info);
    if (sameface==TAO_TRUE) {
      break;
    } else {
      //      info = TT1->WhichBetween(XL,X,XU); CHKERRQ(info);
      TT3=TT2;
      TT2=TT1;
      TT1=TT3;
    }

  }
  
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TronCheckOptimalFace"
static int TronCheckOptimalFace(TaoVec *X,TaoVec *XL,TaoVec*XU,TaoVec *PG,TaoVec*W,
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
    info = tron->TT->whichNonNegative(W); CHKERRQ(info);
    info = tron->TT->GetSize(&n); CHKERRQ(info);
    if (n==0) *optimal = TAO_TRUE;
  */
  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoDefaultMonitor_TRON" 
int TaoDefaultMonitor_TRON(TAO_SOLVER tao,void *dummy)
{
  int info;
  TaoInt its,nfree,nbind;
  double fct,gnorm;
  TAO_TRON *tron;

  TaoFunctionBegin;
  info = TaoGetSolutionStatus(tao,&its,&fct,&gnorm,0,0,0);CHKERRQ(info);
  info = TaoGetSolverContext(tao,"tao_tron",(void**)&tron); CHKERRQ(info);
  if (tron){
    nfree=tron->n_free;
    nbind=tron->n_bind;
    info=TaoPrintInt(tao,"iter = %d,",its); CHKERRQ(info);
    info=TaoPrintDouble(tao," Function value: %g,",fct); CHKERRQ(info);
    info=TaoPrintDouble(tao,"  Residual: %g \n",gnorm);CHKERRQ(info);
    
    info=TaoPrintInt(tao," free vars = %d,",nfree); CHKERRQ(info);
    info=TaoPrintInt(tao," binding vars = %d\n",nbind); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

int TaoSetMaxGPIts(TAO_SOLVER tao, int its){
  int info;
  TAO_TRON  *tron;

  TaoFunctionBegin;

  info = TaoGetSolverContext(tao,"tao_tron",(void**)&tron); CHKERRQ(info);
  if (tron){
    tron->maxgpits     = its;
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetDualVariables_TRON" 
static int TaoGetDualVariables_TRON(TAO_SOLVER tao, TaoVec* DXL, TaoVec* DXU, void *solver){

  TAO_TRON *tron = (TAO_TRON *) solver;
  TaoVec  *G=tron->G,*GP=tron->Work;
  TaoVec  *X,*XL,*XU;
  int       info;

  TaoFunctionBegin;
  info = TaoGetSolution(tao,&X); CHKERRQ(info);
  info = TaoGetVariableBounds(tao,&XL,&XU); CHKERRQ(info);
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
#define __FUNCT__ "TaoCreate_TRON"
int TaoCreate_TRON(TAO_SOLVER tao)
{
  TAO_TRON *tron;
  int      info;

  TaoFunctionBegin;

  info = TaoNew(TAO_TRON,&tron); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_TRON)); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_TRON,(void*)tron); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_TRON,TaoSetDown_TRON); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetOptions_TRON); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_TRON); CHKERRQ(info);
  info=TaoSetTaoDualVariablesRoutine(tao,TaoGetDualVariables_TRON); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao, 50); CHKERRQ(info);
  info = TaoSetTolerances(tao, 1e-10, 1e-10, 0, 0); CHKERRQ(info);

  info = TaoSetTrustRegionRadius(tao, 1.0); CHKERRQ(info);
  info = TaoSetTrustRegionTolerance(tao, 1.0e-12); CHKERRQ(info);

  /* Initialize pointers and variables */
  tron->n            = 0;
  tron->delta        = -1.0;
  tron->maxgpits     = 3;
  tron->pg_ftol      = 0.001;

  tron->eta1         = 1.0e-4;
  tron->eta2         = 0.25;
  tron->eta3         = 0.50;
  tron->eta4         = 0.90;

  tron->sigma1       = 0.5;
  tron->sigma2       = 2.0;
  tron->sigma3       = 4.0;

  tron->gp_iterates  = 0; /* Cumulative number */
  tron->cgits        = 0; /* Current iteration */
  tron->total_gp_its = 0;
  tron->cg_iterates  = 0;
  tron->total_cgits  = 0;
 
  tron->n_bind       = 0;
  tron->n_free       = 0;
  tron->n_upper      = 0;
  tron->n_lower      = 0;

  tron->DX=0;
  tron->DXFree=0;
  tron->R=0;
  tron->X_New=0;
  tron->G_New=0;
  tron->Work=0;
  tron->Free_Local=0;
  tron->TT=0;
  tron->Hsub=0;

  info = TaoCreateMoreThuenteBoundLineSearch(tao, 0 , 0); CHKERRQ(info);
  TaoFunctionReturn(0);
}
EXTERN_C_END
