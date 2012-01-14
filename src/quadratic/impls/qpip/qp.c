#include "qp.h"

static int QPIPSolveLinearSystem(KSP,Mat,Mat,Vec,int,double,Vec,int *);
static int VecAdjust(Vec, Vec);


/* --------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNC__  
#define __FUNC__ "TaoCreate_QPIP"
int TaoCreate_QPIP(TAO_SOLVER tao)
{
  TAO_QPIP *qp;
  KSP       ksp;
  PC        pc;
  int       info;

  PetscFunctionBegin;

  tao->setup              = TaoSetUp_QPIP;
  tao->solve              = TaoSolve_QPIP;
  tao->destroy            = TaoDestroy_QPIP;
  tao->converged          = TaoConverged_Default;
  tao->printhelp          = TaoPrintHelp_QPIP;
  tao->view               = TaoView_QPIP;
  tao->setfromoptions     = TaoSetOptions_QPIP;
  tao->defaultmonitor     = TaoDefaultMonitor_QPIP;
  tao->nwork              = 0;
  tao->fatol               = 1e-8;
  tao->frtol               = 0.0;
  tao->gatol               = 0.0;
  tao->grtol               = 0.0;
  tao->catol               = 1e-8;
  tao->crtol               = 0.0;
  tao->max_its            = 100;
  tao->max_funcs          = 500;
  tao->linear_its         = 0;
  tao->iter               = 0;

  qp                      = PetscNew(TAO_QPIP);CHKPTRQ(qp);
  PLogObjectMemory(tao,sizeof(TAO_QPIP));
  tao->data               = (void *) qp;

  /* Initialize pointers and variables */
  qp->n              = 0;
  qp->m              = 0;
  qp->ksp_tol       = 0.1;
  qp->dobj           = 0.0;
  qp->pobj           = 1.0;
  qp->gap            = 10.0;
  qp->rgap           = 1.0;
  qp->mu             = 1.0;
  qp->sigma          = 1.0;
  qp->dinfeas        = 1.0;
  qp->predcorr       = 1;
  qp->psteplength    = 0.0;
  qp->dsteplength    = 0.0;
  /* Create the linear solver.  By default, use conjugate gradient
     method, block jacobi precondtioner */

  info = KSPCreate(PETSC_COMM_WORLD,&tao->ksp); CHKERRQ(info);  

  info = KSPSetType(ksp,KSPBICG); CHKERRQ(info);
  info = KSPSetInitialGuessNonzero(ksp);CHKERRQ(info);
  info = KSPGetPC(tao->ksp,&pc); CHKERRQ(info);
  info = PCSetType(pc,PCJACOBI); CHKERRQ(info);

  info = KSPSetFromOptions(tao->ksp); CHKERRQ(info);

  info = PetscObjectComposeFunctionDynamic((PetscObject)tao,
                                    "TaoGetKSP_C",
                                    "TaoGetKSP_Default",
                                    (void*)TaoGetKSP_Default);
  CHKERRQ(info);
  /*
  info = PetscObjectComposeFunctionDynamic((PetscObject)tao,
                                    "TaoSetUpperBounds_C",
                                    "TaoSetUpperBounds_Default",
                                    (void*)TaoSetUpperBounds_Default);
  CHKERRQ(info);
  info = PetscObjectComposeFunctionDynamic((PetscObject)tao,
                                    "TaoSetLowerBounds_C",
                                    "TaoSetLowerBounds_Default",
                                    (void*)TaoSetLowerBounds_Default);
  CHKERRQ(info);
  */
  /*
  info = PetscObjectComposeFunctionDynamic((PetscObject)tao,
                                    "TaoGetDualVariables_C",
                                    "TaoGetDualVariables_QPIP",
                                    (void*)TaoGetDualVariables_QPIP);
  CHKERRQ(info);
  */
  PetscFunctionReturn(0);
}
EXTERN_C_END


#undef __FUNC__  
#define __FUNC__ "TaoSetUp_QPIP"
static int TaoSetUp_QPIP(TAO_SOLVER solver)
{
  TAO_QPIP *qp =(TAO_QPIP*)solver->data;
  int n,info;
  Vec XY=qp->XY,XL=qp->XL,XU=qp->XU,D=qp->D,F=qp->F;
  IS ISTemp;
  IS MatOrder[2];

  PLogEventDeactivate(MAT_SetValues);

  info=TaoCheckBounds(solver);CHKERRQ(info);
  info=TaoCheckFG(solver);CHKERRQ(info);
  info=TaoGetVariableBounds(solver,&qp->XL,&qp->XU);CHKERRQ(info);

  /*
  info=TaoGetSolution(solver,&qp->X);CHKERRQ(info);
  info=TaoGetGradient(solver,&qp->Grad);CHKERRQ(info);
  */

  info=VecGetSize(qp->X,&qp->n);CHKERRQ(info);



  info=MatCreateSubspaceProjection(XY, qp->X,qp->ISx,&qp->ProjX); 
  CHKERRQ(info);
  info=MatCreateSubspaceProjection(XY, qp->XU,qp->ISxu,&qp->ProjXU); 
  CHKERRQ(info);
  info=MatCreateSubspaceProjection(XY, qp->XL,qp->ISxl,&qp->ProjXL); 
  CHKERRQ(info);
  info=MatCreateSubspaceProjection(XY, qp->D,qp->ISd,&qp->ProjCL); 
  CHKERRQ(info);
  info=MatCreateSubspaceProjection(XY, qp->F,qp->ISf,&qp->ProjCU); 
  CHKERRQ(info);
  info=MatCreateSubspaceProjection(XY, qp->B,qp->ISEQ,&qp->ProjEQ); 
  CHKERRQ(info);

  info = VecWhichEqual(XY,XY,&ISTemp);CHKERRQ(info);
  info = MatCreateSubspaceProjection(XY, XY,ISTemp,&qp->HPC); CHKERRQ(info);
  info = ISDestroy(ISTemp);CHKERRQ(info);

  /* Allocate some arrays */
  info=VecDuplicate(XY,&qp->DXY); CHKERRQ(info);

  info=VecDuplicate(XL,&qp->G); CHKERRQ(info);
  info=VecDuplicate(XL,&qp->DG); CHKERRQ(info);
  info=VecDuplicate(XL,&qp->Z); CHKERRQ(info);
  info=VecDuplicate(XL,&qp->DZ); CHKERRQ(info);
  info=VecDuplicate(XL,&qp->GZwork); CHKERRQ(info);
  info=VecDuplicate(XL,&qp->R3); CHKERRQ(info);

  info=VecDuplicate(D,&qp->W); CHKERRQ(info);
  info=VecDuplicate(D,&qp->DW); CHKERRQ(info);
  info=VecDuplicate(D,&qp->V); CHKERRQ(info);
  info=VecDuplicate(D,&qp->DV); CHKERRQ(info);
  info=VecDuplicate(D,&qp->WVwork); CHKERRQ(info);
  info=VecDuplicate(D,&qp->R4); CHKERRQ(info);

  info=VecDuplicate(XU,&qp->T); CHKERRQ(info);
  info=VecDuplicate(XU,&qp->DT); CHKERRQ(info);
  info=VecDuplicate(XU,&qp->S); CHKERRQ(info);
  info=VecDuplicate(XU,&qp->DS); CHKERRQ(info);
  info=VecDuplicate(XU,&qp->TSwork); CHKERRQ(info);
  info=VecDuplicate(XU,&qp->R5); CHKERRQ(info);

  info=VecDuplicate(F,&qp->P); CHKERRQ(info);
  info=VecDuplicate(F,&qp->DP); CHKERRQ(info);
  info=VecDuplicate(F,&qp->Q); CHKERRQ(info);
  info=VecDuplicate(F,&qp->DQ); CHKERRQ(info);
  info=VecDuplicate(F,&qp->PQwork); CHKERRQ(info);
  info=VecDuplicate(F,&qp->R6); CHKERRQ(info);

  info=VecDuplicate(XY,&qp->R12); CHKERRQ(info);
  info=VecDuplicate(XY,&qp->RHS); CHKERRQ(info);
  info=VecDuplicate(XY,&qp->RHS2); CHKERRQ(info);
  info=VecDuplicate(XY,&qp->Work); CHKERRQ(info);
  info=VecDuplicate(XY,&qp->DiagAxpy); CHKERRQ(info);
  info=VecDuplicate(XY,&qp->HDiag); CHKERRQ(info);
  info=VecDuplicate(XY,&qp->DScale); CHKERRQ(info);
  info=VecDuplicate(XY,&qp->DScaleH1); CHKERRQ(info);
  info=VecDuplicate(XY,&qp->DScaleH2); CHKERRQ(info);

  info=VecDuplicate(qp->B,&qp->Ywork); CHKERRQ(info);

  qp->m=0;
  info=VecGetSize(qp->XL,&n);  CHKERRQ(info); qp->m+=n;
  info=VecGetSize(qp->XU,&n);  CHKERRQ(info); qp->m+=n;
  info=VecGetSize(qp->D,&n);  CHKERRQ(info); qp->m+=n;
  info=VecGetSize(qp->F,&n);  CHKERRQ(info); qp->m+=n;
 
  info=MatMult(qp->ProjX,qp->C0,qp->C);CHKERRQ(info);
 {
  int n1,n2,n3,n4,n5,n6;
  info=VecGetSize(qp->X,&n1); CHKERRQ(info);
  info=VecGetSize(qp->XL,&n2); CHKERRQ(info);
  info=VecGetSize(qp->XU,&n3); CHKERRQ(info);
  info=VecGetSize(qp->D,&n4); CHKERRQ(info);
  info=VecGetSize(qp->F,&n5); CHKERRQ(info);
  info=VecGetSize(qp->B,&n6); CHKERRQ(info);

  PetscPrintf(PETSC_COMM_WORLD,"\n Variables: %d,  Lower Bounds: %d",n1,n2);
  PetscPrintf(PETSC_COMM_WORLD," Upper Bounds: %d \n",n3);
  PetscPrintf(PETSC_COMM_WORLD," Inequalities : %d,%d, Equality: %d\n",n4,n5,n6);
  }
 /*
 info=MatCreateCmp(solver->comm,2,MatOrder, Vec * X,
		 2, MatOrder, Vec *Y, &qp->H); CHKERRQ(info);
 */
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "QPIPSetInitialPoint"
static int  QPIPSetInitialPoint(TAO_SOLVER tao)
{
  TAO_QPIP *qp = (TAO_QPIP*)tao->data;
  int info;
  double minusone=-1.0,zero=0.0,one=1.0,onehalf=0.5,p01=0.001;
  double gap1,gap2,gap3,gap4,fff,d1,vnorm,mu;

  info = MatMultTranspose(qp->ProjX,qp->X,qp->XY); CHKERRQ(info);

  /* Compute function, Gradient R=Hx+b, and Hessian */

  info = MatMult(qp->H,qp->XY,qp->R12); CHKERRQ(info);
  info = VecDot(qp->C,qp->X,&d1); CHKERRQ(info);
  info = VecDot(qp->R12,qp->XY,&fff); CHKERRQ(info);
  qp->pobj = fff/2.0 + d1 + qp->c;
  tao->fc=qp->pobj;

  /* Initialize Primal Vectors */
  info = MatMult(qp->ProjXL,qp->XY,qp->G);CHKERRQ(info);
  info = VecAXPY(&minusone,qp->XL,qp->G);CHKERRQ(info);
  info = VecSet(&p01,qp->GZwork); CHKERRQ(info);
  info = VecPointwiseMax(qp->G,qp->GZwork,qp->G); CHKERRQ(info);
  
  info = MatMult(qp->ProjCL,qp->R12,qp->W);CHKERRQ(info);
  info = VecScale(&minusone,qp->W);CHKERRQ(info);
  info = VecAXPY(&minusone,qp->D,qp->W);CHKERRQ(info);
  info = VecSet(&p01,qp->WVwork); CHKERRQ(info);
  info = VecPointwiseMax(qp->W,qp->WVwork,qp->W); CHKERRQ(info);
  
  info = MatMult(qp->ProjXU,qp->XY,qp->T);CHKERRQ(info);
  info = VecAXPY(&minusone,qp->XU,qp->T);CHKERRQ(info);
  info = VecScale(&minusone,qp->T);CHKERRQ(info);
  info = VecSet(&p01,qp->TSwork); CHKERRQ(info);
  info = VecPointwiseMax(qp->T,qp->TSwork,qp->T); CHKERRQ(info);
  
  info = MatMult(qp->ProjCU,qp->R12,qp->P);CHKERRQ(info);
  info = VecAXPY(&one,qp->F,qp->P);CHKERRQ(info);
  info = VecSet(&p01,qp->PQwork); CHKERRQ(info);
  info = VecPointwiseMax(qp->P,qp->PQwork,qp->P); CHKERRQ(info);

  /* Initialize Dual Variable Vectors */
  info = VecCopy(qp->G,qp->Z); CHKERRQ(info);
  info = VecReciprocal(qp->Z); CHKERRQ(info);
  info = VecCopy(qp->T,qp->S); CHKERRQ(info);
  info = VecReciprocal(qp->S); CHKERRQ(info);
  info = VecCopy(qp->W,qp->V); CHKERRQ(info);
  info = VecReciprocal(qp->V); CHKERRQ(info);
  info = VecCopy(qp->P,qp->Q); CHKERRQ(info);
  info = VecReciprocal(qp->Q); CHKERRQ(info);

  info = MatMult(qp->H,qp->Work,qp->RHS); CHKERRQ(info);
  info = VecSelectiveSign(qp->ISEQ, qp->DScale);CHKERRQ(info);
  info = VecAXPY(&minusone,qp->DScale,qp->RHS); CHKERRQ(info);
  info = VecAbs(qp->RHS); CHKERRQ(info);
  info = VecSet(&p01,qp->Work); CHKERRQ(info);
  info = VecPointwiseMax(qp->RHS,qp->Work,qp->RHS); CHKERRQ(info);

  info = VecSet(&zero,qp->Work); CHKERRQ(info);
  info = VecSelectiveCopy(qp->Work, qp->ISEQ, qp->R12); CHKERRQ(info);
  info = VecPointwiseDivide(qp->R12,qp->RHS,qp->RHS); CHKERRQ(info);
  info = VecNorm(qp->RHS,NORM_1,&vnorm); CHKERRQ(info);
  mu = TaoMin(1.0,(vnorm+10.0)/qp->m);

  info = VecScale(&mu,qp->S); CHKERRQ(info);
  info = VecScale(&mu,qp->Z); CHKERRQ(info);
  info = VecScale(&mu,qp->V); CHKERRQ(info);
  info = VecScale(&mu,qp->Q); CHKERRQ(info);

  /* Compute the duality gap */
  info = VecDot(qp->G,qp->Z, &gap1); CHKERRQ(info);
  info = VecDot(qp->T,qp->S, &gap2); CHKERRQ(info);
  info = VecDot(qp->W,qp->V, &gap3); CHKERRQ(info);
  info = VecDot(qp->P,qp->Q, &gap4); CHKERRQ(info);
  qp->gap = (gap1+gap2+gap3+gap4);
  
  info = QPIPComputeResidual(qp); CHKERRQ(info);
  info = TAOComputeNormFromCentralPath_QPIP(tao,&d1); CHKERRQ(info);
  /* Compute the dual infeasibility; R = Hx+c+s-z  */
  qp->rnorm=(qp->dinfeas+qp->pinfeas)/(qp->m+qp->n);
  qp->dobj = qp->pobj - qp->gap;
  tao->fc=qp->pobj;
  if (qp->m>0) qp->mu=qp->gap/(qp->m); else qp->mu=0.0;
  qp->rgap=qp->gap/( fabs(qp->dobj) + fabs(qp->pobj) + 1.0 );
  tao->norm=qp->gap + qp->dinfeas;
  tao->cnorm0=qp->pinfeas;
  
  while ( qp->rnorm>3*qp->mu && qp->mu >0) {

    info=TaoDefaultMonitor_QPIP(tao,TAO_NULL);CHKERRQ(info);
    info=VecShift(&qp->mu,qp->G); CHKERRQ(info);
    info=VecShift(&qp->mu,qp->Z); CHKERRQ(info);
    info=VecShift(&qp->mu,qp->S); CHKERRQ(info);
    info=VecShift(&qp->mu,qp->T); CHKERRQ(info);
    info=VecShift(&qp->mu,qp->W); CHKERRQ(info);
    info=VecShift(&qp->mu,qp->V); CHKERRQ(info);
    info=VecShift(&qp->mu,qp->P); CHKERRQ(info);
    info=VecShift(&qp->mu,qp->Q); CHKERRQ(info);

    info = QPIPComputeResidual(qp); CHKERRQ(info);
    
    /* Compute the duality gap */
    info = VecDot(qp->G,qp->Z, &gap1); CHKERRQ(info);
    info = VecDot(qp->T,qp->S, &gap2); CHKERRQ(info);
    info = VecDot(qp->W,qp->V, &gap3); CHKERRQ(info);
    info = VecDot(qp->P,qp->Q, &gap4); CHKERRQ(info);
    
    qp->gap = (gap1+gap2+gap3+gap4);
    qp->dobj = qp->pobj - qp->gap;
    tao->fc=qp->pobj;
    if (qp->m>0) qp->mu=qp->gap/(qp->m); else qp->mu=0.0;
    qp->rgap=qp->gap/( fabs(qp->dobj) + fabs(qp->pobj) + 1.0 );
    tao->norm=qp->gap + qp->dinfeas;
    tao->cnorm0=qp->pinfeas;
    info = TAOComputeNormFromCentralPath_QPIP(tao,&d1); CHKERRQ(info);
    qp->rnorm=(qp->dinfeas+qp->pinfeas)/(qp->m+qp->n);
  }

  PetscFunctionReturn(0);
}


#undef __FUNC__  
#define __FUNC__ "TaoDestroy_QPIP"
static int TaoDestroy_QPIP(TAO_SOLVER tao)
{
  TAO_QPIP *qp = (TAO_QPIP*)tao->data;
  KSP ksp;
  int info;
  /* Free allocated memory in GPCG structure */
  
  PetscFunctionBegin;

  if (tao->setupcalled) {
    info=VecDestroy(qp->DXY); CHKERRQ(info); 
    info=VecDestroy(qp->R12); CHKERRQ(info);
    info=VecDestroy(qp->HDiag); CHKERRQ(info);
    info=VecDestroy(qp->Work); CHKERRQ(info);
    info=VecDestroy(qp->XY); CHKERRQ(info);
    info=VecDestroy(qp->DiagAxpy); CHKERRQ(info);
    info=VecDestroy(qp->DScale); CHKERRQ(info);
    info=VecDestroy(qp->RHS); CHKERRQ(info);
    info=VecDestroy(qp->RHS2); CHKERRQ(info);
    
    info=VecDestroy(qp->G); CHKERRQ(info);
    info=VecDestroy(qp->DG); CHKERRQ(info);
    info=VecDestroy(qp->Z); CHKERRQ(info);
    info=VecDestroy(qp->DZ); CHKERRQ(info);
    info=VecDestroy(qp->GZwork); CHKERRQ(info);
    info=VecDestroy(qp->R3); CHKERRQ(info);
    
    info=VecDestroy(qp->V); CHKERRQ(info);
    info=VecDestroy(qp->W); CHKERRQ(info);
    info=VecDestroy(qp->DW); CHKERRQ(info);
    info=VecDestroy(qp->DV); CHKERRQ(info);
    info=VecDestroy(qp->WVwork); CHKERRQ(info);
    info=VecDestroy(qp->R4); CHKERRQ(info);
    
    info=VecDestroy(qp->S); CHKERRQ(info);
    info=VecDestroy(qp->DS); CHKERRQ(info);
    info=VecDestroy(qp->T); CHKERRQ(info);
    info=VecDestroy(qp->DT); CHKERRQ(info);
    info=VecDestroy(qp->TSwork); CHKERRQ(info);
    info=VecDestroy(qp->R5); CHKERRQ(info);
    
    info=VecDestroy(qp->P); CHKERRQ(info);
    info=VecDestroy(qp->DP); CHKERRQ(info);
    info=VecDestroy(qp->Q); CHKERRQ(info);
    info=VecDestroy(qp->DQ); CHKERRQ(info);
    info=VecDestroy(qp->PQwork); CHKERRQ(info);
    info=VecDestroy(qp->R6); CHKERRQ(info);
    
    info=VecDestroy(qp->Ywork); CHKERRQ(info);
  }

  info=TaoGetKSP(tao,&ksp);CHKERRQ(info);
  info=KSPDestroy(ksp);CHKERRQ(info);
  PetscFree(qp);

  PetscFunctionReturn(0);
}


#undef __FUNC__  
#define __FUNC__ "TaoSolve_QPIP"
static int TaoSolve_QPIP(TAO_SOLVER solver, int *its)
{
  TAO_QPIP *qp = (TAO_QPIP*)solver->data;
  int info,cgits,maxcgiters=qp->n;
  double minusone=-1.0,zero=0.0,one=1.0;
  double ksptol,sigma, sigmamu;
  double dstep,pstep,step,d1,d2,gap[4]; 
  TaoTruth flg;
  TaoTerminateReason reason=TAO_CONTINUE_ITERATING;

  PetscFunctionBegin;

  info=TaoCheckBounds(solver);CHKERRQ(info);
  info=TaoCheckFG(solver);CHKERRQ(info);

  info = OptionsGetLogical(TAO_NULL,"-predcorr",&qp->predcorr, &flg);
  CHKERRQ(info);

  if (solver->userstart==TAO_FALSE){
  }

  info = MatGetDiagonal(qp->H,qp->HDiag); CHKERRQ(info);
  info = MatMaxPerRow(qp->H,qp->DScaleH1); CHKERRQ(info);
  info = VecAbs(qp->DScaleH1);CHKERRQ(info);
  info = VecSet(&one,qp->Work);CHKERRQ(info);
  info = VecPointwiseMax(qp->DScaleH1,qp->Work,qp->DScaleH1);CHKERRQ(info);

  info = QPIPSetInitialPoint(solver); CHKERRQ(info);
 
  solver->norm=qp->dinfeas;
  solver->norm0=qp->dinfeas;
  solver->cnorm=qp->pinfeas;
  solver->cnorm0=qp->pinfeas;

  for ( qp->iter = 0; qp->iter < solver->max_its+2; (qp->iter)++ ){

    info = TAOComputeNormFromCentralPath_QPIP(solver,&d1); CHKERRQ(info);
    info = TaoDefaultMonitor_QPIP(solver,TAO_NULL);CHKERRQ(info);
    TaoMonitor(solver);
    TaoLogConvHistory(solver,solver->norm,solver->iter);
    
    /* Check Stopping Condition      */
    info = TaoCheckConvergence(solver);CHKERRQ(info);
    info = TaoGetTerminationReason(solver,&reason); CHKERRQ(info);
    if (reason != TAO_CONTINUE_ITERATING) break;
    
    /* 
       Dual Infeasibility Direction should already be in the right
       hand side from computing the residuals 
    */
    
    if (qp->iter < 0 && (qp->rnorm>5*qp->mu || d1>sqrt(qp->m)*qp->mu) ) {
      sigma=1.0;sigmamu=qp->mu;
      printf("Toward Central Path \n");      
    } else {
      sigma=0.0;sigmamu=0;
    }
    info = VecSet(&sigmamu,qp->DZ); CHKERRQ(info);
    info = VecSet(&sigmamu,qp->DS); CHKERRQ(info);
    info = VecSet(&sigmamu,qp->DV); CHKERRQ(info);
    info = VecSet(&sigmamu,qp->DQ); CHKERRQ(info);
    
    if (sigmamu !=0){
      info = VecPointwiseDivide(qp->DZ,qp->G,qp->DZ); CHKERRQ(info);
      info = VecPointwiseDivide(qp->DS,qp->T,qp->DS); CHKERRQ(info);
      info = VecPointwiseDivide(qp->DV,qp->W,qp->DV); CHKERRQ(info);
      info = VecPointwiseDivide(qp->DQ,qp->P,qp->DQ); CHKERRQ(info);

      info = MatMultTranspose(qp->ProjXL,qp->DZ,qp->RHS2);CHKERRQ(info);
      info = MatMultTranspose(qp->ProjXU,qp->DS,qp->Work);CHKERRQ(info);
      info = VecAXPY(&minusone,qp->Work,qp->RHS2); CHKERRQ(info);
      info = MatMultTranspose(qp->ProjCL,qp->DV,qp->Work);CHKERRQ(info);
      info = VecAXPY(&one,qp->Work,qp->RHS2); CHKERRQ(info);
      info = MatMultTranspose(qp->ProjCU,qp->DQ,qp->Work);CHKERRQ(info);
      info = VecAXPY(&minusone,qp->Work,qp->RHS2); CHKERRQ(info);

    } else {
      info = VecSet(&zero,qp->RHS2); CHKERRQ(info);
    }

    info = VecPointwiseDivide(qp->Z,qp->G,qp->GZwork); CHKERRQ(info);
    info = MatMultTranspose(qp->ProjXL,qp->GZwork,qp->DiagAxpy);CHKERRQ(info);
    info = VecPointwiseMult(qp->GZwork,qp->R3,qp->GZwork); CHKERRQ(info);
    info = MatMultTranspose(qp->ProjXL,qp->GZwork,qp->Work);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->Work,qp->RHS); CHKERRQ(info);

    info = VecPointwiseDivide(qp->S,qp->T,qp->TSwork); CHKERRQ(info);
    info = MatMultTranspose(qp->ProjXU,qp->TSwork,qp->Work);CHKERRQ(info);
    info = VecAXPY(&one,qp->Work,qp->DiagAxpy); CHKERRQ(info);
    info = VecPointwiseMult(qp->TSwork,qp->R5,qp->TSwork); CHKERRQ(info);
    info = MatMultTranspose(qp->ProjXU,qp->TSwork,qp->Work);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->Work,qp->RHS); CHKERRQ(info);

    info = VecPointwiseDivide(qp->V,qp->W,qp->WVwork); CHKERRQ(info);
    info = MatMultTranspose(qp->ProjCL,qp->WVwork,qp->DScale);CHKERRQ(info);
    info = VecPointwiseMult(qp->WVwork,qp->R4,qp->WVwork); CHKERRQ(info);
    info = MatMultTranspose(qp->ProjCL,qp->WVwork,qp->Work);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->Work,qp->RHS); CHKERRQ(info);

    info = VecPointwiseDivide(qp->Q,qp->P,qp->PQwork); CHKERRQ(info);
    info = MatMultTranspose(qp->ProjCU,qp->PQwork,qp->Work);CHKERRQ(info);
    info = VecAXPY(&one,qp->Work,qp->DScale); CHKERRQ(info);
    info = VecPointwiseMult(qp->PQwork,qp->R6,qp->PQwork); CHKERRQ(info);
    info = MatMultTranspose(qp->ProjCU,qp->PQwork,qp->Work);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->Work,qp->RHS); CHKERRQ(info);

    info = VecAXPY(&one,qp->RHS,qp->RHS2); CHKERRQ(info);
    info = VecAdjust(qp->DScale,qp->DiagAxpy); CHKERRQ(info);

    info = VecPointwiseMult(qp->RHS2,qp->DScale,qp->RHS2); CHKERRQ(info);
    info = VecPointwiseMult(qp->DiagAxpy,qp->DScale,qp->DiagAxpy); CHKERRQ(info);

    /*  Determine the solving tolerance */
    ksptol = qp->mu/10;
    ksptol = TaoMin(ksptol,0.001);

    /* Approximately solve the linear system */
    info = MatDiagonalShift(qp->H,qp->DiagAxpy); CHKERRQ(info);

    info = MatGetDiagonal(qp->H,qp->Work); CHKERRQ(info);
    info = VecAbs(qp->Work); CHKERRQ(info);
    info = VecSet(&one,qp->DScaleH1); CHKERRQ(info);
    info = VecPointwiseMax(qp->DScaleH1,qp->Work,qp->Work);CHKERRQ(info);

    info = MatSetDiagonal(qp->HPC,qp->Work); CHKERRQ(info);

    /*    info = MatView(qp->H,VIEWER_DRAW_WORLD); CHKERRQ(info); */
    info = QPIPSolveLinearSystem(solver->ksp,qp->H,qp->H,qp->RHS2,
				 maxcgiters,ksptol,
				 qp->DXY,&cgits); CHKERRQ(info);

    info = MatSetDiagonal(qp->H,qp->HDiag); CHKERRQ(info);
    info = QPComputeStepDirection(qp); CHKERRQ(info);
    info = QPStepLength(qp); CHKERRQ(info);

    pstep = qp->psteplength; dstep = qp->dsteplength;
    
    info = VecDot(qp->DG,qp->DZ, gap); CHKERRQ(info);
    info = VecDot(qp->DT,qp->DS, gap+1); CHKERRQ(info);
    info = VecDot(qp->DW,qp->DV, gap+2); CHKERRQ(info);
    info = VecDot(qp->DP,qp->DQ, gap+3); CHKERRQ(info);
    
    sigmamu= ( pstep*pstep*(gap[0]+gap[1]+gap[2]+gap[3]) +
	       (1 - pstep + pstep*sigma)*qp->gap  )/qp->m;
    
    if (qp->predcorr && (step < 0.9 || sigma==1.0) ){
      if (sigma==1.0){
	sigmamu=1.0;
      } else if (sigmamu < qp->mu){ 
	sigmamu=sigmamu/qp->mu;
	sigmamu=sigmamu*sigmamu*sigmamu;
      } else {sigmamu = 1.0;}
      sigmamu = sigmamu*qp->mu;
      /* Compute Corrector Step */
      info = VecPointwiseMult(qp->DG,qp->DZ,qp->DZ); CHKERRQ(info);
      info = VecScale(&minusone,qp->DZ); CHKERRQ(info);
      info = VecShift(&sigmamu,qp->DZ); CHKERRQ(info);
      info = VecPointwiseDivide(qp->DZ,qp->G,qp->DZ); CHKERRQ(info);
      
      info = VecPointwiseMult(qp->DW,qp->DV,qp->DV); CHKERRQ(info);
      info = VecScale(&minusone,qp->DV); CHKERRQ(info);
      info = VecShift(&sigmamu,qp->DV); CHKERRQ(info);
      info = VecPointwiseDivide(qp->DV,qp->W,qp->DV); CHKERRQ(info);
      
      info = VecPointwiseMult(qp->DT,qp->DS,qp->DS); CHKERRQ(info);
      info = VecScale(&minusone,qp->DS); CHKERRQ(info);
      info = VecShift(&sigmamu,qp->DS); CHKERRQ(info);
      info = VecPointwiseDivide(qp->DS,qp->T,qp->DS); CHKERRQ(info);

      info = VecPointwiseMult(qp->DP,qp->DQ,qp->DQ); CHKERRQ(info);
      info = VecScale(&minusone,qp->DQ); CHKERRQ(info);
      info = VecShift(&sigmamu,qp->DQ); CHKERRQ(info);
      info = VecPointwiseDivide(qp->DQ,qp->P,qp->DQ); CHKERRQ(info);

      info = MatMultTranspose(qp->ProjXL,qp->DZ,qp->RHS2);CHKERRQ(info);
      info = MatMultTranspose(qp->ProjXU,qp->DS,qp->Work);CHKERRQ(info);
      info = VecAXPY(&minusone,qp->Work,qp->RHS2); CHKERRQ(info);
      info = MatMultTranspose(qp->ProjCL,qp->DV,qp->Work);CHKERRQ(info);
      info = VecAXPY(&one,qp->Work,qp->RHS2); CHKERRQ(info);
      info = MatMultTranspose(qp->ProjCU,qp->DQ,qp->Work);CHKERRQ(info);
      info = VecAXPY(&minusone,qp->Work,qp->RHS2); CHKERRQ(info);

      info = VecAXPY(&one,qp->RHS,qp->RHS2); CHKERRQ(info);
      info = VecPointwiseMult(qp->RHS2,qp->DScale,qp->RHS2); CHKERRQ(info);
      
      /*  Determine the solving tolerance */
      
      ksptol = sigmamu;
      ksptol = TaoMin(ksptol,0.0001);

      info = MatDiagonalShift(qp->H, qp->DiagAxpy); CHKERRQ(info);
      
      info = QPIPSolveLinearSystem(solver->ksp,qp->H,qp->H,qp->RHS2,
				   maxcgiters,ksptol,
				   qp->DXY,&cgits); CHKERRQ(info);
      
      info = VecScale(&minusone,qp->DiagAxpy); CHKERRQ(info);
      info = MatDiagonalShift(qp->H,qp->DiagAxpy); CHKERRQ(info);
      info = VecScale(&minusone,qp->DiagAxpy); CHKERRQ(info);
      info = QPComputeStepDirection(qp); CHKERRQ(info);
      info = QPStepLength(qp);  CHKERRQ(info);

    }  /* End Corrector step */

    if (1==0){
      pstep = 1.0; dstep = 1.0;
      
      info = VecAXPY(&dstep, qp->DZ, qp->Z);  CHKERRQ(info);
      info = VecAXPY(&dstep, qp->DS, qp->S);  CHKERRQ(info);
      info = VecAXPY(&dstep, qp->DV, qp->V);  CHKERRQ(info);
      info = VecAXPY(&dstep, qp->DQ, qp->Q);  CHKERRQ(info);
      info = VecAXPY(&dstep, qp->DXY, qp->XY);  CHKERRQ(info);
      info = VecAXPY(&pstep, qp->DG, qp->G);  CHKERRQ(info);
      info = VecAXPY(&pstep, qp->DT, qp->T);  CHKERRQ(info);
      info = VecAXPY(&pstep, qp->DW, qp->W);  CHKERRQ(info);
      info = VecAXPY(&pstep, qp->DP, qp->P);  CHKERRQ(info);
      
      /* Compute Residuals */
      info = QPIPComputeResidual(qp); CHKERRQ(info);
      
      pstep = -1.0; dstep = -1.0;
      
      info = VecAXPY(&dstep, qp->DZ, qp->Z);  CHKERRQ(info);
      info = VecAXPY(&dstep, qp->DS, qp->S);  CHKERRQ(info);
      info = VecAXPY(&dstep, qp->DV, qp->V);  CHKERRQ(info);
      info = VecAXPY(&dstep, qp->DQ, qp->Q);  CHKERRQ(info);
      info = VecAXPY(&dstep, qp->DXY, qp->XY);  CHKERRQ(info);
      info = VecAXPY(&pstep, qp->DG, qp->G);  CHKERRQ(info);
      info = VecAXPY(&pstep, qp->DT, qp->T);  CHKERRQ(info);
      info = VecAXPY(&pstep, qp->DW, qp->W);  CHKERRQ(info);
      info = VecAXPY(&pstep, qp->DP, qp->P);  CHKERRQ(info);
    }
    
    pstep = qp->psteplength; dstep = qp->dsteplength;

    info = VecAXPY(&dstep, qp->DZ, qp->Z);  CHKERRQ(info);
    info = VecAXPY(&dstep, qp->DS, qp->S);  CHKERRQ(info);
    info = VecAXPY(&dstep, qp->DV, qp->V);  CHKERRQ(info);
    info = VecAXPY(&dstep, qp->DQ, qp->Q);  CHKERRQ(info);
    info = VecAXPY(&dstep, qp->DXY, qp->XY);  CHKERRQ(info);
    info = VecAXPY(&pstep, qp->DG, qp->G);  CHKERRQ(info);
    info = VecAXPY(&pstep, qp->DT, qp->T);  CHKERRQ(info);
    info = VecAXPY(&pstep, qp->DW, qp->W);  CHKERRQ(info);
    info = VecAXPY(&pstep, qp->DP, qp->P);  CHKERRQ(info);
    
    /* Compute Residuals */
    info = QPIPComputeResidual(qp); CHKERRQ(info);

    /* Update the vectors */
    info = MatMult(qp->ProjX,qp->XY,qp->X); CHKERRQ(info);
    info = MatMult(qp->ProjX,qp->R12,qp->Grad); CHKERRQ(info);
    /*    info = MatMult(qp->ProjX,qp->DXY,qp->DX); CHKERRQ(info); */
    /*    info = MatMult(qp->ProjY,qp->XY,qp->Y); CHKERRQ(info); */

    /* Evaluate quadratic function */
    info = MatMultTranspose(qp->ProjX,qp->X,qp->DiagAxpy); CHKERRQ(info);
    info = MatMult(qp->H,qp->DiagAxpy,qp->Work); CHKERRQ(info);
    info = VecDot(qp->Work,qp->DiagAxpy,&d1); CHKERRQ(info);
    info = VecDot(qp->X,qp->C,&d2); CHKERRQ(info);
    qp->pobj=d1/2.0 + d2 + qp->c;

    /* Compute the duality gap */
    info = VecDot(qp->W,qp->V, gap); CHKERRQ(info);
    info = VecDot(qp->P,qp->Q, gap+1); CHKERRQ(info);
    info = VecDot(qp->G,qp->Z, gap+2); CHKERRQ(info);
    info = VecDot(qp->T,qp->S, gap+3); CHKERRQ(info);

    solver->norm=qp->dinfeas+qp->gap;
    solver->cnorm=qp->pinfeas;
    qp->gap = (gap[0]+gap[1]+gap[2]+gap[3]);
    qp->dobj = qp->pobj - qp->gap;
    if (qp->m>0) qp->mu=qp->gap/(qp->m);
    qp->rgap=qp->gap/( fabs(qp->dobj) + fabs(qp->pobj) + 1.0 );
    solver->fc=qp->pobj;

    info = TaoIncrementIterate(solver);CHKERRQ(info);

  }  /* END MAIN LOOP  */
  
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "QPComputeStepDirection"
static int QPComputeStepDirection(TAO_QPIP *qp){

  int info;
  double minusone = -1.0, one = 1.0;

  PetscFunctionBegin;
  info = MatMult(qp->H,qp->DXY,qp->Work); CHKERRQ(info);
    
  /* Calculate DG */
  info = MatMult(qp->ProjXL,qp->DXY,qp->DG);CHKERRQ(info);
  info = VecAXPY(&one,qp->R3,qp->DG); CHKERRQ(info);

  /* Calculate DZ */  
  info = VecAXPY(&minusone,qp->Z,qp->DZ);  CHKERRQ(info);  
  info = VecPointwiseDivide(qp->DG,qp->G,qp->GZwork); CHKERRQ(info);
  info = VecPointwiseMult(qp->GZwork,qp->Z,qp->GZwork); CHKERRQ(info);
  info = VecAXPY(&minusone,qp->GZwork,qp->DZ); CHKERRQ(info);

  /* Calculate DT */
  info = MatMult(qp->ProjXU,qp->DXY,qp->DT);CHKERRQ(info);
  info = VecScale(&minusone,qp->DT); CHKERRQ(info);
  info = VecAXPY(&minusone,qp->R5,qp->DT); CHKERRQ(info);

  /* Calculate DS */
  info = VecAXPY(&minusone, qp->S, qp->DS);  CHKERRQ(info);
  info = VecPointwiseDivide(qp->DT,qp->T,qp->TSwork); CHKERRQ(info);
  info = VecPointwiseMult(qp->TSwork,qp->S,qp->TSwork); CHKERRQ(info);
  info = VecAXPY(&minusone, qp->TSwork, qp->DS);  CHKERRQ(info);

  /* Calculate DW */
  info = MatMult(qp->ProjCL,qp->Work,qp->DW);CHKERRQ(info);
  info = VecScale(&minusone,qp->DW); CHKERRQ(info);
  info = VecAXPY(&one,qp->R4,qp->DW); CHKERRQ(info);

  /* Calculate DV */
  info = VecAXPY(&minusone,qp->V,qp->DV);  CHKERRQ(info);
  info = VecPointwiseDivide(qp->DW,qp->W,qp->WVwork); CHKERRQ(info);
  info = VecPointwiseMult(qp->WVwork,qp->V,qp->WVwork); CHKERRQ(info);
  info = VecAXPY(&minusone,qp->WVwork,qp->DV); CHKERRQ(info);
  
  /* Calculate DP */
  info = MatMult(qp->ProjCU,qp->Work,qp->DP);CHKERRQ(info);
  info = VecAXPY(&minusone,qp->R6,qp->DP); CHKERRQ(info);

  /* Calculate DQ */  
  info = VecAXPY(&minusone, qp->Q, qp->DQ);  CHKERRQ(info);
  info = VecPointwiseDivide(qp->DP,qp->P,qp->PQwork); CHKERRQ(info);
  info = VecPointwiseMult(qp->PQwork,qp->Q,qp->PQwork); CHKERRQ(info);
  info = VecAXPY(&minusone, qp->PQwork, qp->DQ);  CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "QPIPComputeResidual"
static int QPIPComputeResidual(TAO_QPIP *qp){

  int info;
  double minusone = -1.0, one = 1.0;
  double R12norm,R3norm,R4norm,R5norm,R6norm,R7norm;
  double dtmp = 1.0 - qp->psteplength;

  PetscFunctionBegin;
  info = MatMult(qp->H,qp->XY,qp->RHS); CHKERRQ(info);
  info = MatMultTranspose(qp->ProjX,qp->C,qp->Work); CHKERRQ(info);
  info = VecAXPY(&one,qp->Work,qp->RHS); CHKERRQ(info);
  info = VecScale(&minusone,qp->RHS); CHKERRQ(info);

  info = MatMultTranspose(qp->ProjEQ,qp->B,qp->Work);CHKERRQ(info);
  info = VecAXPY(&minusone,qp->Work,qp->RHS);CHKERRQ(info);

  if (1==0){

    info = VecScale(&dtmp,qp->R3); CHKERRQ(info);
    info = VecScale(&dtmp,qp->R4); CHKERRQ(info);
    info = VecScale(&dtmp,qp->R5); CHKERRQ(info);
    info = VecScale(&dtmp,qp->R6); CHKERRQ(info);

  } else {

    info = MatMult(qp->ProjXL,qp->XY,qp->R3);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->XL,qp->R3);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->G,qp->R3);CHKERRQ(info);
    
    info = MatMult(qp->ProjCL,qp->RHS,qp->R4);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->D,qp->R4);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->W,qp->R4);CHKERRQ(info);
    
    info = MatMult(qp->ProjXU,qp->XY,qp->R5);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->XU,qp->R5);CHKERRQ(info);
    info = VecAXPY(&one,qp->T,qp->R5);CHKERRQ(info);
    
    info = MatMult(qp->ProjCU,qp->RHS,qp->R6);CHKERRQ(info);
    info = VecAXPY(&minusone,qp->F,qp->R6);CHKERRQ(info);
    info = VecAXPY(&one,qp->P,qp->R6);CHKERRQ(info);
    
  }


  /* Compute R7 */
  info = MatMult(qp->ProjEQ,qp->RHS,qp->Ywork);CHKERRQ(info);

  info = VecSelectiveCopy(qp->XY,qp->ISd,qp->RHS);CHKERRQ(info);
  info = VecSelectiveSign(qp->ISd,qp->RHS);CHKERRQ(info);
  info = VecSelectiveCopy(qp->XY,qp->ISf,qp->RHS);CHKERRQ(info);
  info = VecSelectiveSign(qp->ISf,qp->RHS);CHKERRQ(info);

  info = MatMultTranspose(qp->ProjXU,qp->S,qp->R12);CHKERRQ(info);
  info = MatMultTranspose(qp->ProjXL,qp->Z,qp->Work);CHKERRQ(info);
  info = VecAXPY(&minusone,qp->Work,qp->R12);CHKERRQ(info);
  info = MatMultTranspose(qp->ProjCU,qp->Q,qp->Work);CHKERRQ(info);
  info = VecAXPY(&one,qp->Work,qp->R12);CHKERRQ(info);
  info = MatMultTranspose(qp->ProjCL,qp->V,qp->Work);CHKERRQ(info);
  info = VecAXPY(&minusone,qp->Work,qp->R12);CHKERRQ(info);

  info = VecAXPY(&minusone,qp->RHS,qp->R12);CHKERRQ(info);

  info = VecNorm(qp->R12,NORM_2, &R12norm); CHKERRQ(info);
  info = VecNorm(qp->R3,NORM_2, &R3norm); CHKERRQ(info);
  info = VecNorm(qp->R4,NORM_2, &R4norm); CHKERRQ(info);
  info = VecNorm(qp->R5,NORM_2, &R5norm); CHKERRQ(info);
  info = VecNorm(qp->R6,NORM_2, &R6norm); CHKERRQ(info);
  info = VecNorm(qp->Ywork,NORM_2, &R7norm); CHKERRQ(info);

  qp->rnorm=R12norm;
  qp->dinfeas=sqrt(R12norm*R12norm-R7norm*R7norm);

  qp->pinfeas=TaoMax(R3norm,R4norm);
  qp->pinfeas=TaoMax(R4norm,qp->pinfeas);
  qp->pinfeas=TaoMax(R5norm,qp->pinfeas);
  qp->pinfeas=TaoMax(R6norm,qp->pinfeas);
  qp->pinfeas=TaoMax(R7norm,qp->pinfeas);

  if (0==1){  
    printf("R12 = %4.2e, ",R12norm);
    printf("R3 = %4.2e, ",R3norm);
    printf("R4 = %4.2e, ",R4norm);
    printf("R5 = %4.2e, ",R5norm);
    printf("R6 = %4.2e, ",R6norm);
    printf("R7 = %4.2e,\n ",R7norm);
    
    info = MatMult(qp->ProjXL,qp->R12,qp->GZwork);CHKERRQ(info);
    info = MatMult(qp->ProjXU,qp->R12,qp->TSwork);CHKERRQ(info);
    info = MatMult(qp->ProjCL,qp->R12,qp->WVwork);CHKERRQ(info);
    info = MatMult(qp->ProjCU,qp->R12,qp->PQwork);CHKERRQ(info);
    
    info = VecNorm(qp->GZwork,NORM_2, &R3norm); CHKERRQ(info);
    info = VecNorm(qp->TSwork,NORM_2, &R4norm); CHKERRQ(info);
    info = VecNorm(qp->WVwork,NORM_2, &R5norm); CHKERRQ(info);
    info = VecNorm(qp->PQwork,NORM_2, &R6norm); CHKERRQ(info);
    printf("R3 = %4.2e, ",R3norm);
    printf("R4 = %4.2e, ",R4norm);
    printf("R5 = %4.2e, ",R5norm);
    printf("R6 = %4.2e\n",R6norm);
  }

  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "QPStepLength"
static int QPStepLength(TAO_QPIP *qp){

  double tstep1,tstep2,tstep3,tstep4,tstep=1.0;
  double np5 = -0.5;
  double scl=0.99;
  int info;

  PetscFunctionBegin;
  /* Compute stepsize to the boundary */
  if (qp->mu>1.0) scl=0.9;
  else if (qp->mu>0.01) scl=0.95;
  else scl=0.99;

  info = VecPointwiseDivide(qp->DG,qp->G,qp->GZwork); CHKERRQ(info);
  info = VecBoundAbove(&np5,qp->GZwork); CHKERRQ(info);

  info = VecPointwiseDivide(qp->DT,qp->T,qp->TSwork); CHKERRQ(info);
  info = VecBoundAbove(&np5,qp->TSwork); CHKERRQ(info);

  info = VecPointwiseDivide(qp->DW,qp->W,qp->WVwork); CHKERRQ(info);
  info = VecBoundAbove(&np5,qp->WVwork); CHKERRQ(info);

  info = VecPointwiseDivide(qp->DP,qp->P,qp->PQwork); CHKERRQ(info);
  info = VecBoundAbove(&np5,qp->PQwork); CHKERRQ(info);

  info = VecNorm(qp->GZwork,NORM_INFINITY, &tstep1); CHKERRQ(info);
  info = VecNorm(qp->TSwork,NORM_INFINITY, &tstep2); CHKERRQ(info);
  info = VecNorm(qp->WVwork,NORM_INFINITY, &tstep3); CHKERRQ(info);
  info = VecNorm(qp->PQwork,NORM_INFINITY, &tstep4); CHKERRQ(info);

  tstep = TaoMax(tstep1,tstep2);
  tstep = TaoMax(tstep,tstep3);
  tstep = TaoMax(tstep,tstep4);
  tstep = TaoMin(1.0,scl/tstep);
  qp->psteplength = tstep;

  info = VecPointwiseDivide(qp->DS,qp->S,qp->TSwork); CHKERRQ(info);
  info = VecBoundAbove(&np5,qp->TSwork); CHKERRQ(info);

  info = VecPointwiseDivide(qp->DZ,qp->Z,qp->GZwork); CHKERRQ(info);
  info = VecBoundAbove(&np5,qp->GZwork); CHKERRQ(info);

  info = VecPointwiseDivide(qp->DQ,qp->Q,qp->PQwork); CHKERRQ(info);
  info = VecBoundAbove(&np5,qp->PQwork); CHKERRQ(info);

  info = VecPointwiseDivide(qp->DV,qp->V,qp->WVwork); CHKERRQ(info);
  info = VecBoundAbove(&np5,qp->WVwork); CHKERRQ(info);

  info = VecNorm(qp->GZwork,NORM_INFINITY, &tstep1); CHKERRQ(info);
  info = VecNorm(qp->TSwork,NORM_INFINITY, &tstep2); CHKERRQ(info);
  info = VecNorm(qp->WVwork,NORM_INFINITY, &tstep3); CHKERRQ(info);
  info = VecNorm(qp->PQwork,NORM_INFINITY, &tstep4); CHKERRQ(info);

  tstep = TaoMax(tstep1,tstep2);
  tstep = TaoMax(tstep,tstep3);
  tstep = TaoMax(tstep,tstep4);
  tstep = TaoMin(1.0,scl/tstep);
  qp->dsteplength = tstep;

  qp->psteplength = TaoMin(qp->psteplength,qp->dsteplength);
  qp->dsteplength = qp->psteplength;

  PetscFunctionReturn(0);
}



#undef __FUNC__  
#define __FUNC__ "QPIPSolveLinearSystem"
static int QPIPSolveLinearSystem(KSP ksp,Mat P,Mat Q,Vec R,int maxits,
				 double tol,Vec DXFree,int *its){

  int info;
  KSP ksp;
  double zero=0.0;
  double abstol=10e-12, dtol = 1.0e30;

  PetscFunctionBegin;
  /* Set solver operators and options */
  info = KSPSetOperators(ksp,Q,P,DIFFERENT_NONZERO_PATTERN);

  CHKERRQ(info);

  info = KSPSetTolerances(ksp,zero,tol,dtol,maxits/3+2); CHKERRQ(info);

  info = KSPSolve(ksp,R,DXFree,its); CHKERRQ(info);
  if (*its<0){
    info = KSPSolve(ksp,R,DXFree,&its2); CHKERRQ(info);
    if (its2<0) its=its+its2; else *its=its2- *its;
  }
  if (*its<0){
    info = KSPSolve(ksp,R,DXFree,&its2); CHKERRQ(info);
    if (*its2<0) *its=*its+its2; else *its=its2-*its;
  }
  PetscFunctionReturn(0);
}


#undef __FUNC__  
#define __FUNC__ "TAOComputeNormFromCentralPath_QPIP"
int TAOComputeNormFromCentralPath_QPIP(TAO_SOLVER solver, double *norm)
{
  TAO_QPIP *qp = (TAO_QPIP*)solver->data;
  int       info;
  double    gap[4], nmu;
  
  PetscFunctionBegin;
  
  info = VecPointwiseMult(qp->G,qp->Z,qp->GZwork); CHKERRQ(info);  
  info = VecPointwiseMult(qp->T,qp->S,qp->TSwork); CHKERRQ(info);
  info = VecPointwiseMult(qp->W,qp->V,qp->WVwork); CHKERRQ(info);
  info = VecPointwiseMult(qp->P,qp->Q,qp->PQwork); CHKERRQ(info);

  nmu=-qp->mu;

  info = VecShift(&nmu,qp->GZwork); CHKERRQ(info);
  info = VecShift(&nmu,qp->TSwork); CHKERRQ(info);
  info = VecShift(&nmu,qp->WVwork); CHKERRQ(info);
  info = VecShift(&nmu,qp->PQwork); CHKERRQ(info);

  info = VecDotBegin(qp->GZwork,qp->GZwork,&gap[0]); CHKERRQ(info);
  info = VecDotBegin(qp->TSwork,qp->TSwork,&gap[1]); CHKERRQ(info);
  info = VecDotBegin(qp->WVwork,qp->WVwork,&gap[2]); CHKERRQ(info);
  info = VecDotBegin(qp->PQwork,qp->PQwork,&gap[3]); CHKERRQ(info);
  info = VecDotEnd(qp->GZwork,qp->GZwork,&gap[0]); CHKERRQ(info);
  info = VecDotEnd(qp->TSwork,qp->TSwork,&gap[1]); CHKERRQ(info);
  info = VecDotEnd(qp->WVwork,qp->WVwork,&gap[2]); CHKERRQ(info);
  info = VecDotEnd(qp->PQwork,qp->PQwork,&gap[3]); CHKERRQ(info);
  
  qp->pathnorm=sqrt( (gap[0]+gap[1]+gap[2]+gap[3]) );
  *norm=qp->pathnorm;

  PetscFunctionReturn(0);
}



int VecSelectiveSign(IS Select, Vec V){

  int i,n,in,*ip,low,high,info;
  double *v;

  info = VecGetOwnershipRange(V, &low, &high); CHKERRQ(info);

  info = VecGetLocalSize(V,&n); CHKERRQ(info);
  info = ISGetSize(Select, &in); CHKERRQ(info);

  info = VecGetArray(V,&v); CHKERRQ(info);

  info = ISGetIndices(Select, &ip); CHKERRQ(info);

  for (i=0; i<in; i++){
    if (ip[i]>= low && ip[i]<high)
      v[ip[i]-low] *= -1.0 ;
    else{
      PetscPrintf(PETSC_COMM_WORLD,
		  " %d of %d index: %2d   low: %d  high: %d \n",
		  i,in,ip[i],low,high); 
      SETERRQ(1,0,"IS index out of range");
    }
  }

  info = ISRestoreIndices(Select,&ip); CHKERRQ(info);

  info = VecRestoreArray(V,&v); CHKERRQ(info);

  return 0;
}

/*
int VecISAXPY(double * alpha, IS XIS, Vec X, IS YIS, Vec Y){

  int i,info,in1,in2,xlow,xhigh,ylow,yhigh,*xis,*yis;
  double *x,*y;


  info = ISGetSize(XIS, &in1); CHKERRQ(info);
  info = ISGetSize(YIS, &in2); CHKERRQ(info);
  if ( in1 != in2 )
    SETERRQ(1,0,"Index sets must be identically loaded over processors");

  info = VecGetOwnershipRange(X, &xlow, &xhigh); CHKERRQ(info);
  info = VecGetOwnershipRange(Y, &ylow, &yhigh); CHKERRQ(info);

  info = ISGetIndices(XIS, &xis); CHKERRQ(info);
  if (XIS != YIS){
    info = ISGetIndices(YIS, &yis); CHKERRQ(info);
  } else {
    yis=xis;
  }

  info = VecGetArray(X,&x); CHKERRQ(info);
  if (X != Y){
    info = VecGetArray(Y,&y); CHKERRQ(info);
  } else {
    y=x;
  }

  for (i=0; i<in1; i++){
    if (yis[i]>=ylow && yis[i]<yhigh && xis[i]>=xlow && xis[i]<xhigh){
      y[yis[i]-ylow] = y[yis[i]-ylow] + (*alpha)*x[xis[i]-xlow];
    } else {
      SETERRQ(1,0,"IS index out of range");
    }      
  }

  info = ISRestoreIndices(XIS, &xis); CHKERRQ(info);
  if (XIS != YIS){
    info = ISRestoreIndices(YIS, &yis); CHKERRQ(info);
  }

  info = VecRestoreArray(X,&x); CHKERRQ(info);
  if (X != Y){
    info = VecRestoreArray(Y,&y); CHKERRQ(info);
  }

  PLogFlops(2*in1);

  return 0;
}
*/

int VecBoundAbove(double * b, Vec V){

  int i,n,info;
  double *vptr;

  info = VecGetLocalSize(V,&n); CHKERRQ(info);

  info = VecGetArray(V,&vptr); CHKERRQ(info);
  for (i=0; i<n; i++){
    vptr[i] = TaoMin(vptr[i],*b);
  }
  info = VecRestoreArray(V,&vptr); CHKERRQ(info);

  return 0;
}



#undef __FUNC__  
#define __FUNC__ "TaoDefaultMonitor_QPIP"
static int TaoDefaultMonitor_QPIP(TAO_SOLVER tao,void *y)
{
  TAO_QPIP *qp = (TAO_QPIP*)tao->data;
  int       info;
  double    dtemp;

  PetscFunctionBegin;
  info=TAOComputeNormFromCentralPath_QPIP(tao,&dtemp);
  CHKERRQ(info);

  info = PetscPrintf(tao->comm,"iter = %d, Fcn value: %g, Infeas: %g, ",
		     tao->iter,qp->pobj, qp->pinfeas);CHKERRQ(info);
  info = PetscPrintf(tao->comm," Dinfeas: %g, Gap:%g,Nbhood:%g Mu:%g\n",
		     qp->dinfeas,qp->gap,dtemp,qp->mu);CHKERRQ(info);

  info = PetscPrintf(tao->comm,"step: %g\n",qp->psteplength);CHKERRQ(info);
  PetscFunctionReturn(0);
}


static int TaoSetOptions_QPIP(TAO_SOLVER solver){

  TAO_QPIP *qp = (TAO_QPIP*)solver->data;
  int       info;
  TaoTruth  flg;

  PetscFunctionBegin;
  info = OptionsGetLogical(TAO_NULL,"-predcorr",&qp->predcorr, &flg);
  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "TaoView_QPIP"
static int TaoView_QPIP(TAO_SOLVER tao,Viewer v){
  int info;
  KSP taoksp;

  PetscFunctionBegin;
  info = TaoGetKSP(tao,&taoksp);CHKERRQ(info);
  info = KSPView(taoksp,v);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "TaoPrintHelp_QPIP"
static int TaoPrintHelp_QPIP(TAO_SOLVER t,char * c)
{ 
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}


int SetBQPIPOperators(Mat Q, Vec X, IS ISx, double c, Vec C0,
                      Vec XL, IS ISxl, Vec XU, IS ISxu, Vec B, IS ISEQ,
                      Vec D, IS ISd, Vec F, IS ISf, TAO_SOLVER tao){

  /* Set pointers to Data */
  TAO_QPIP *qp = (TAO_QPIP*)tao->data;
  int info;

  qp->XL=XL;
  qp->XU=XU;
  qp->H=Q;
  qp->C0=C0;
  qp->B=B;
  qp->c=-c;
  qp->D=D;
  qp->F=F;

  info=VecDuplicate(C0,&qp->XY);CHKERRQ(info);
  info=VecDuplicate(X,&qp->C);CHKERRQ(info);
  info=VecDuplicate(X,&qp->X);CHKERRQ(info);
  info=VecDuplicate(X,&qp->Grad);CHKERRQ(info);
  info=VecDuplicate(X,&qp->DX);CHKERRQ(info);

  qp->ISx = ISx;
  qp->ISxl = ISxl;
  qp->ISxu = ISxu;
  qp->ISEQ = ISEQ;
  qp->ISd = ISd;
  qp->ISf = ISf;

  return 0;
}



static int VecAdjust(Vec Scale, Vec Diag){

  int i,n,info;
  double *sptr,*dptr;

  info = VecGetArray(Scale,&sptr);CHKERRQ(info);
  info = VecGetArray(Diag,&dptr);CHKERRQ(info);
  info = VecGetLocalSize(Scale,&n);CHKERRQ(info);
  for (i=0;i<n;i++){
    if (sptr[i]>0){
      sptr[i]=-1.0/sptr[i];
      dptr[i]=1.0;
    } else {
      sptr[i]=1.0;
    }
  }

  info = VecRestoreArray(Scale,&sptr);CHKERRQ(info);
  info = VecRestoreArray(Diag,&dptr);CHKERRQ(info);

  return 0;
}




#include "src/vec/vecimpl.h"    /*I "tao_solver.h" I*/
#undef __FUNC__  
#define __FUNC__ "VecCreateProductSpace"
int VecCreateProductSpace(Vec A, Vec B, Vec *J, IS*ISA, IS* ISB){

  int info,ione=1;
  int n1,n2,nl1,nl2;
  int la,lb,ha,hb;
  Scalar zero=0.0;
  VecType type_name;

  PetscFunctionBegin;
  PetscCheckSameComm(A,B);

  info = VecGetSize(A,&n1);CHKERRQ(info);
  info = VecGetSize(B,&n2);CHKERRQ(info);
  info = VecGetLocalSize(A,&nl1);CHKERRQ(info);
  info = VecGetLocalSize(B,&nl2);CHKERRQ(info);
  info = VecGetOwnershipRange(A,&la,&ha);CHKERRQ(info);
  info = VecGetOwnershipRange(B,&lb,&hb);CHKERRQ(info);
  if (J){
    info = VecCreate(A->comm,nl1+nl2,n1+n2,J);CHKERRQ(info);
    info = VecGetType(A,&type_name);CHKERRQ(info);
    info = VecSetType(*J,type_name);CHKERRQ(info);
    info = VecSet(&zero,*J);CHKERRQ(info);
  }
  if (ISA){
    info = ISCreateStride(A->comm,ha-la,la+lb,ione,ISA);CHKERRQ(info);
  }
  if (ISB){
    info = ISCreateStride(A->comm,hb-lb,ha+lb,ione,ISB);CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}






