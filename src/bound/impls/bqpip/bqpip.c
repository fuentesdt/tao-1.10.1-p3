/*$Id$*/

#include "bqpip.h"
#include "taolinearsolver.h"
static int TaoSetDown_BQPIP(TAO_SOLVER, void*);

#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_BQPIP"
static int TaoSetUp_BQPIP(TAO_SOLVER tao, void *solver)
{
  TAO_BQPIP *qp =(TAO_BQPIP*)solver;
  TaoInt       n, info;
  TaoLinearSolver *ksp;
  TaoFunctionBegin;

  /* Set pointers to Data */
  info = TaoGetHessian(tao,&qp->H);
  info = TaoGetSolution(tao,&qp->XY);CHKERRQ(info);
  info = qp->XY->GetDimension(&qp->n); CHKERRQ(info);

  /* Allocate some arrays */
  info = qp->XY->Clone(&qp->Work); CHKERRQ(info);
  info = qp->XY->Clone(&qp->DXY); CHKERRQ(info);
  info = qp->XY->Clone(&qp->HDiag); CHKERRQ(info);
  info = qp->XY->Clone(&qp->DiagAxpy); CHKERRQ(info);
  info = qp->XY->Clone(&qp->RHS); CHKERRQ(info);
  info = qp->XY->Clone(&qp->RHS2); CHKERRQ(info);
  info = qp->XY->Clone(&qp->C0); CHKERRQ(info);
  info = qp->XY->Clone(&qp->R12); CHKERRQ(info);

  info = qp->XY->Clone(&qp->XL); CHKERRQ(info);
  info = qp->XL->Clone(&qp->G); CHKERRQ(info);
  info = qp->XL->Clone(&qp->DG); CHKERRQ(info);
  info = qp->XL->Clone(&qp->Z); CHKERRQ(info);
  info = qp->XL->Clone(&qp->DZ); CHKERRQ(info);
  info = qp->XL->Clone(&qp->GZwork); CHKERRQ(info);
  info = qp->XL->Clone(&qp->R3); CHKERRQ(info);

  info = qp->XY->Clone(&qp->XU); CHKERRQ(info);
  info = qp->XU->Clone(&qp->T); CHKERRQ(info);
  info = qp->XU->Clone(&qp->DT); CHKERRQ(info);
  info = qp->XU->Clone(&qp->S); CHKERRQ(info);
  info = qp->XU->Clone(&qp->DS); CHKERRQ(info);
  info = qp->XU->Clone(&qp->TSwork); CHKERRQ(info);
  info = qp->XU->Clone(&qp->R5); CHKERRQ(info);

  info = TaoSetLagrangianGradientVector(tao,qp->R12);CHKERRQ(info);
  info = TaoSetVariableBounds(tao,qp->XL,qp->XU);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,qp->DXY);CHKERRQ(info);

  /* Register the Events */
  info = qp->XY->GetDimension(&qp->n); CHKERRQ(info);

  qp->m=0;
  info=qp->G->GetDimension(&n); CHKERRQ(info); qp->m+=n;
  info=qp->T->GetDimension(&n); CHKERRQ(info); qp->m+=n;

  info = TaoCreateLinearSolver(tao,qp->H,300,&ksp); CHKERRQ(info);
  info = ksp->SetTolerances(1e-14,1e-30,1e30,qp->n); CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPIPSetInitialPoint"
static int  QPIPSetInitialPoint(TAO_SOLVER tao, TAO_BQPIP *qp)
{
  int       info;
  double    two=2.0,p01=1;
  double    gap1,gap2,fff,mu;

  TaoFunctionBegin;
  /* Compute function, Gradient R=Hx+b, and Hessian */
  info = qp->XY->Median(qp->XL,qp->XY,qp->XU); CHKERRQ(info);
  info = qp->H->Multiply(qp->XY,qp->R12); CHKERRQ(info);

  info = qp->Work->Waxpby(0.5,qp->R12,1.0,qp->C0);CHKERRQ(info);
  info = qp->R12->Axpy(1.0,qp->C0);CHKERRQ(info);
  info = qp->Work->Dot(qp->XY,&fff);CHKERRQ(info);
  qp->pobj = fff + qp->c;

  /* Initialize Primal Vectors */

  info = qp->T->Waxpby(1.0,qp->XU,-1.0,qp->XY);CHKERRQ(info);
  info = qp->G->Waxpby(1.0,qp->XY,-1.0,qp->XL);CHKERRQ(info);

  info = qp->GZwork->SetToConstant(p01);CHKERRQ(info);
  info = qp->TSwork->SetToConstant(p01);CHKERRQ(info);

  info = qp->G->PointwiseMaximum(qp->G,qp->GZwork); CHKERRQ(info);
  info = qp->T->PointwiseMaximum(qp->T,qp->TSwork); CHKERRQ(info);

  /* Initialize Dual Variable Vectors */

  info = qp->Z->CopyFrom(qp->G); CHKERRQ(info);
  info = qp->Z->Reciprocal(); CHKERRQ(info);

  info = qp->S->CopyFrom(qp->T); CHKERRQ(info);
  info = qp->S->Reciprocal(); CHKERRQ(info);

  info = qp->H->Multiply(qp->Work,qp->RHS); CHKERRQ(info);
  info = qp->RHS->AbsoluteValue(); CHKERRQ(info);
  info = qp->Work->SetToConstant(p01);CHKERRQ(info);
  info = qp->RHS->PointwiseMaximum(qp->RHS,qp->Work); CHKERRQ(info);

  info = qp->RHS->PointwiseDivide(qp->R12,qp->RHS); CHKERRQ(info);
  info = qp->RHS->Norm1(&gap1); CHKERRQ(info);
  mu = TaoMin(10.0,(gap1+10.0)/qp->m);

  info = qp->S->Scale(mu); CHKERRQ(info);
  info = qp->Z->Scale(mu); CHKERRQ(info);

  info = qp->TSwork->SetToConstant(p01); CHKERRQ(info);
  info = qp->GZwork->SetToConstant(p01); CHKERRQ(info);
  info = qp->S->PointwiseMaximum(qp->S,qp->TSwork); CHKERRQ(info);
  info = qp->Z->PointwiseMaximum(qp->Z,qp->GZwork); CHKERRQ(info);

  qp->mu=0;qp->dinfeas=1.0;qp->pinfeas=1.0;
  while ( (qp->dinfeas+qp->pinfeas)/(qp->m+qp->n) >= qp->mu ){

    info=qp->G->Scale(two); CHKERRQ(info);
    info=qp->Z->Scale(two); CHKERRQ(info);
    info=qp->S->Scale(two); CHKERRQ(info);
    info=qp->T->Scale(two); CHKERRQ(info);

    info = QPIPComputeResidual(qp); CHKERRQ(info);
    
    info=qp->R3->Waxpby(1.0,qp->XY,-1.0,qp->G);CHKERRQ(info);
    info=qp->R3->Axpy(-1.0,qp->XL);CHKERRQ(info);

    info=qp->R5->Waxpby(1.0,qp->XY,1.0,qp->T);CHKERRQ(info);
    info=qp->R5->Axpy(-1.0,qp->XU);CHKERRQ(info);
    
    info=qp->R3->NormInfinity(&gap1);CHKERRQ(info);
    info=qp->R5->NormInfinity(&gap2);CHKERRQ(info);
    qp->pinfeas=TaoMax(gap1,gap2);
    
    /* Compute the duality gap */
    info=qp->G->Dot(qp->Z,&gap1);CHKERRQ(info);
    info=qp->T->Dot(qp->S,&gap2);CHKERRQ(info);
    
    qp->gap = (gap1+gap2);
    qp->dobj = qp->pobj - qp->gap;
    if (qp->m>0) qp->mu=qp->gap/(qp->m); else qp->mu=0.0;
    qp->rgap=qp->gap/( TaoAbsScalar(qp->dobj) + TaoAbsScalar(qp->pobj) + 1.0 );
  }

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_BQPIP"
static int TaoSetDown_BQPIP(TAO_SOLVER tao, void*solver)
{
  TAO_BQPIP *qp = (TAO_BQPIP*)solver;
  int info;

  /* Free allocated memory in GPCG structure */
  TaoFunctionBegin;

  info=TaoVecDestroy(qp->G);CHKERRQ(info);
  info=TaoVecDestroy(qp->DG);CHKERRQ(info);
  info=TaoVecDestroy(qp->Z);CHKERRQ(info);
  info=TaoVecDestroy(qp->DZ);CHKERRQ(info);
  info=TaoVecDestroy(qp->GZwork);CHKERRQ(info);
  info=TaoVecDestroy(qp->R3);CHKERRQ(info);
  info=TaoVecDestroy(qp->XL);CHKERRQ(info);
  
  info=TaoVecDestroy(qp->S);CHKERRQ(info);
  info=TaoVecDestroy(qp->DS);CHKERRQ(info);
  info=TaoVecDestroy(qp->T);CHKERRQ(info);
  info=TaoVecDestroy(qp->DT);CHKERRQ(info);
  info=TaoVecDestroy(qp->TSwork);CHKERRQ(info);
  info=TaoVecDestroy(qp->R5);CHKERRQ(info);
  info=TaoVecDestroy(qp->XU);CHKERRQ(info);
  
  info=TaoVecDestroy(qp->HDiag);CHKERRQ(info);
  info=TaoVecDestroy(qp->Work);CHKERRQ(info);
  info=TaoVecDestroy(qp->DiagAxpy);CHKERRQ(info);
  info=TaoVecDestroy(qp->RHS);CHKERRQ(info);
  info=TaoVecDestroy(qp->RHS2);CHKERRQ(info);
  info=TaoVecDestroy(qp->DXY);CHKERRQ(info);
  info=TaoVecDestroy(qp->C0);CHKERRQ(info);
  info=TaoVecDestroy(qp->R12);CHKERRQ(info);
  
  info = TaoDestroyLinearSolver(tao);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_BQPIP"
static int TaoSolve_BQPIP(TAO_SOLVER tao, void *solver)
{
  TAO_BQPIP *qp = (TAO_BQPIP*)solver;
  int       info;
  TaoInt    iter=0;
  TaoInt    n;
  double    d1,d2,ksptol,sigma;
  double    sigmamu;
  double    dstep,pstep,step=0;
  double    gap[4];
  TaoTruth  kspsuccess;
  TaoTerminateReason reason;
  
  TaoFunctionBegin;

  info = TaoGetSolution(tao,&qp->XY);CHKERRQ(info);
  info = TaoEvaluateVariableBounds(tao,qp->XL,qp->XU); CHKERRQ(info);
  info = qp->XY->GetDimension(&n);CHKERRQ(info);

  info = TaoComputeFunctionGradient(tao,qp->XY,&qp->c,qp->C0); CHKERRQ(info);
  info = TaoComputeHessian(tao,qp->XY,qp->H); CHKERRQ(info);
  info = qp->H->Multiply(qp->XY,qp->Work); CHKERRQ(info);
  info = qp->XY->Dot(qp->Work,&d1); CHKERRQ(info);
  info = qp->C0->Axpy(-1.0,qp->Work); CHKERRQ(info);
  info = qp->XY->Dot(qp->C0,&d2); CHKERRQ(info);
  qp->c -= (d1/2.0+d2);

  info = qp->H->GetDiagonal(qp->HDiag); CHKERRQ(info);

  info = QPIPSetInitialPoint(tao,qp); CHKERRQ(info);
  info = QPIPComputeResidual(qp); CHKERRQ(info);
  
  /* Enter main loop */
  while (1){

    /* Check Stopping Condition      */
    info=TaoMonitor(tao,iter++,qp->pobj,sqrt(qp->gap+qp->dinfeas),qp->pinfeas,
		    step,&reason); CHKERRQ(info);
    if (reason != TAO_CONTINUE_ITERATING) break;

    /* 
       Dual Infeasibility Direction should already be in the right
       hand side from computing the residuals 
    */

    info = TAOComputeNormFromCentralPath_BQPIP(tao,&d1); CHKERRQ(info);

    if (iter > 0 && (qp->rnorm>5*qp->mu || d1*d1>qp->m*qp->mu*qp->mu) ) {
      sigma=1.0;sigmamu=qp->mu;
      sigma=0.0;sigmamu=0;
    } else {
      sigma=0.0;sigmamu=0;
    }
    info = qp->DZ->SetToConstant(sigmamu); CHKERRQ(info);
    info = qp->DS->SetToConstant(sigmamu); CHKERRQ(info);

    if (sigmamu !=0){
      info = qp->DZ->PointwiseDivide(qp->DZ,qp->G); CHKERRQ(info);
      info = qp->DS->PointwiseDivide(qp->DS,qp->T); CHKERRQ(info);
      info = qp->RHS2->Waxpby(1.0,qp->DZ,1.0,qp->DS); CHKERRQ(info);
    } else {
      info = qp->RHS2->SetToZero(); CHKERRQ(info);
    }


    /* 
       Compute the Primal Infeasiblitiy RHS and the 
       Diagonal Matrix to be added to H and store in Work 
    */
    info = qp->DiagAxpy->PointwiseDivide(qp->Z,qp->G); CHKERRQ(info);
    info = qp->GZwork->PointwiseMultiply(qp->DiagAxpy,qp->R3); CHKERRQ(info);
    info = qp->RHS->Axpy(-1.0,qp->GZwork); CHKERRQ(info);

    info = qp->TSwork->PointwiseDivide(qp->S,qp->T); CHKERRQ(info);
    info = qp->DiagAxpy->Axpy(1.0,qp->TSwork); CHKERRQ(info);
    info = qp->TSwork->PointwiseMultiply(qp->TSwork,qp->R5); CHKERRQ(info);
    info = qp->RHS->Axpy(-1.0,qp->TSwork); CHKERRQ(info);

    info = qp->RHS2->Axpy(1.0,qp->RHS); CHKERRQ(info);

    /*  Determine the solving tolerance */
    ksptol = qp->mu/10.0;
    ksptol = TaoMin(ksptol,0.001);

    info = qp->H->AddDiagonal(qp->DiagAxpy); CHKERRQ(info);

    info = TaoPreLinearSolve(tao,qp->H);CHKERRQ(info);
    info = TaoLinearSolve(tao,qp->H,qp->RHS,qp->DXY,&kspsuccess);CHKERRQ(info);
    

    info = qp->DiagAxpy->Negate(); CHKERRQ(info);
    info = qp->H->AddDiagonal(qp->DiagAxpy); CHKERRQ(info);
    info = qp->DiagAxpy->Negate(); CHKERRQ(info);
    info = QPComputeStepDirection(qp); CHKERRQ(info);
    info = QPStepLength(qp);  CHKERRQ(info);

    /* Calculate New Residual R1 in Work vector */
    info = qp->H->Multiply(qp->DXY,qp->RHS2); CHKERRQ(info);
    info = qp->RHS2->Axpy(1.0,qp->DS); CHKERRQ(info);
    info = qp->RHS2->Axpy(-1.0,qp->DZ); CHKERRQ(info);
    info = qp->RHS2->Aypx(qp->dsteplength,qp->R12); CHKERRQ(info);

    qp->RHS2->Norm2(&qp->dinfeas); CHKERRQ(info);
    qp->DG->Dot(qp->DZ, gap); CHKERRQ(info);
    qp->DT->Dot(qp->DS, gap+1); CHKERRQ(info);
 
    qp->rnorm=(qp->dinfeas+qp->psteplength*qp->pinfeas)/(qp->m+qp->n);
    pstep = qp->psteplength; dstep = qp->dsteplength;    
    step = TaoMin(qp->psteplength,qp->dsteplength);
    sigmamu= ( pstep*pstep*(gap[0]+gap[1]) +
	       (1 - pstep + pstep*sigma)*qp->gap  )/qp->m;

    if (qp->predcorr && step < 0.9){
      if (sigmamu < qp->mu){ 
	sigmamu=sigmamu/qp->mu;
	sigmamu=sigmamu*sigmamu*sigmamu;
      } else {sigmamu = 1.0;}
      sigmamu = sigmamu*qp->mu;
      
      /* Compute Corrector Step */
      info = qp->DZ->PointwiseMultiply(qp->DG,qp->DZ); CHKERRQ(info);
      info = qp->DZ->Negate(); CHKERRQ(info);
      info = qp->DZ->AddConstant(sigmamu); CHKERRQ(info);
      info = qp->DZ->PointwiseDivide(qp->DZ,qp->G); CHKERRQ(info);

      info = qp->DS->PointwiseMultiply(qp->DS,qp->DT); CHKERRQ(info);
      info = qp->DS->Negate(); CHKERRQ(info);
      info = qp->DS->AddConstant(sigmamu); CHKERRQ(info);
      info = qp->DS->PointwiseDivide(qp->DS,qp->T); CHKERRQ(info);
      
      info = qp->RHS2->Waxpby(1.0,qp->DZ,-1.0,qp->DS); CHKERRQ(info);
      info = qp->RHS2->Axpy(1.0,qp->RHS); CHKERRQ(info);

      /* Approximately solve the linear system */
      
      info = qp->H->AddDiagonal(qp->DiagAxpy); CHKERRQ(info);
       
      info = TaoLinearSolve(tao,qp->H,qp->RHS2,qp->DXY,&kspsuccess);CHKERRQ(info);
      
      info = qp->H->SetDiagonal(qp->HDiag); CHKERRQ(info);
      info = QPComputeStepDirection(qp); CHKERRQ(info);
      info = QPStepLength(qp); CHKERRQ(info);

    }  /* End Corrector step */


    /* Take the step */
    pstep = qp->psteplength; dstep = qp->dsteplength;

    info = qp->Z->Axpy(dstep,qp->DZ); CHKERRQ(info);
    info = qp->S->Axpy(dstep,qp->DS); CHKERRQ(info);
    info = qp->XY->Axpy(dstep,qp->DXY); CHKERRQ(info);
    info = qp->G->Axpy(pstep,qp->DG); CHKERRQ(info);
    info = qp->T->Axpy(pstep,qp->DT); CHKERRQ(info);
    
    /* Compute Residuals */
    info = QPIPComputeResidual(qp); CHKERRQ(info);

    /* Evaluate quadratic function */
    info = qp->H->Multiply(qp->XY,qp->Work); CHKERRQ(info);

    info = qp->XY->Dot(qp->Work,&d1); CHKERRQ(info);
    info = qp->XY->Dot(qp->C0,&d2); CHKERRQ(info);
    info = qp->G->Dot(qp->Z,gap); CHKERRQ(info);
    info = qp->T->Dot(qp->S,gap+1); CHKERRQ(info);

    qp->pobj=d1/2.0 + d2+qp->c;
    /* Compute the duality gap */
    qp->gap = (gap[0]+gap[1]);
    qp->dobj = qp->pobj - qp->gap;
    if (qp->m>0) qp->mu=qp->gap/(qp->m);
    qp->rgap=qp->gap/( fabs(qp->dobj) + fabs(qp->pobj) + 1.0 );
  }  /* END MAIN LOOP  */

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPComputeStepDirection"
static int QPComputeStepDirection(TAO_BQPIP *qp)
{
  int info;

  TaoFunctionBegin;

  /* Calculate DG */
  info = qp->DG->Waxpby(1.0,qp->DXY,1.0,qp->R3);CHKERRQ(info);

  /* Calculate DT */
  info = qp->DT->Waxpby(-1.0,qp->DXY,-1.0,qp->R5);CHKERRQ(info);

  /* Calculate DZ */
  
  info = qp->DZ->Axpy(-1.0,qp->Z);CHKERRQ(info);
  info = qp->GZwork->PointwiseDivide(qp->DG,qp->G); CHKERRQ(info);
  info = qp->GZwork->PointwiseMultiply(qp->GZwork,qp->Z); CHKERRQ(info);
  info = qp->DZ->Axpy(-1.0,qp->GZwork);CHKERRQ(info);

  /* Calculate DS */

  info = qp->DS->Axpy(-1.0,qp->S);CHKERRQ(info);
  info = qp->TSwork->PointwiseDivide(qp->DT,qp->T); CHKERRQ(info);
  info = qp->TSwork->PointwiseMultiply(qp->TSwork,qp->S); CHKERRQ(info);
  info = qp->DS->Axpy(-1.0,qp->TSwork);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPIPComputeResidual"
static int QPIPComputeResidual(TAO_BQPIP *qp)
{
  int info;
  double gap1,gap2,dtmp = 1.0 - qp->psteplength;

  TaoFunctionBegin;

  /* Compute R3 and R5 */

  if (1==1){

    info = qp->R3->Scale(dtmp);
    info = qp->R5->Scale(dtmp);
    qp->pinfeas=dtmp*qp->pinfeas;

  } else {

    info = qp->R3->Waxpby(1.0,qp->XY,-1.0,qp->XL);CHKERRQ(info);
    info = qp->R3->Axpy(-1.0,qp->G);CHKERRQ(info);

    info = qp->R5->Waxpby(1.0,qp->XY,-1.0,qp->XU);CHKERRQ(info);
    info = qp->R5->Axpy(1.0,qp->T);CHKERRQ(info);

    info = qp->R3->NormInfinity(&gap1);CHKERRQ(info);
    info = qp->R5->NormInfinity(&gap2);CHKERRQ(info);

    qp->pinfeas=TaoMax(gap1,gap2);

  }

  qp->R12->Waxpby(1.0,qp->S,-1.0,qp->Z);CHKERRQ(info);

  info = qp->H->Multiply(qp->XY,qp->RHS); CHKERRQ(info);
  info = qp->RHS->Negate();CHKERRQ(info);
  info = qp->RHS->Axpy(-1.0,qp->C0);CHKERRQ(info);
  info = qp->R12->Axpy(-1.0,qp->RHS);CHKERRQ(info);

  info = qp->R12->Norm1(&qp->dinfeas); CHKERRQ(info);
  qp->rnorm=(qp->dinfeas+qp->pinfeas)/(qp->m+qp->n);
  
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPStepLength"
static int QPStepLength(TAO_BQPIP *qp)
{
  double tstep1,tstep2,tstep3,tstep4,tstep;
  int info;

  TaoFunctionBegin;
  /* Compute stepsize to the boundary */

  info = qp->G->StepMax(qp->DG,&tstep1); CHKERRQ(info);
  info = qp->T->StepMax(qp->DT,&tstep2); CHKERRQ(info);
  info = qp->S->StepMax(qp->DS,&tstep3); CHKERRQ(info);
  info = qp->Z->StepMax(qp->DZ,&tstep4); CHKERRQ(info);

  tstep = TaoMin(tstep1,tstep2);
  qp->psteplength = TaoMin(0.95*tstep,1.0);

  tstep = TaoMin(tstep3,tstep4);
  qp->dsteplength = TaoMin(0.95*tstep,1.0);

  qp->psteplength = TaoMin(qp->psteplength,qp->dsteplength);
  qp->dsteplength = qp->psteplength;

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoGetDualVariables_BQPIP"
int TaoGetDualVariables_BQPIP(TAO_SOLVER tao,TaoVec* DXL, TaoVec* DXU, void *solver)
{
  TAO_BQPIP *qp = (TAO_BQPIP*)solver;
  int       info;

  TaoFunctionBegin;
  info = DXL->CopyFrom(qp->Z); CHKERRQ(info);
  info = DXU->CopyFrom(qp->S); CHKERRQ(info);
  info = DXU->Negate(); CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TAOComputeNormFromCentralPath_BQPIP"
int TAOComputeNormFromCentralPath_BQPIP(TAO_SOLVER tao, double *norm)
{
  TAO_BQPIP *qp;
  int       info;
  double    gap[2],mu[2], nmu;
  
  TaoFunctionBegin;
  info = TaoGetSolverContext(tao,"tao_bqpip",(void**)&qp); CHKERRQ(info);
  info = qp->GZwork->PointwiseMultiply(qp->G,qp->Z); CHKERRQ(info);
  info = qp->TSwork->PointwiseMultiply(qp->T,qp->S); CHKERRQ(info);
  info = qp->TSwork->Norm1(&mu[0]); CHKERRQ(info);
  info = qp->GZwork->Norm1(&mu[1]); CHKERRQ(info);

  nmu=-(mu[0]+mu[1])/qp->m;

  qp->GZwork->AddConstant(nmu); CHKERRQ(info);
  qp->TSwork->AddConstant(nmu); CHKERRQ(info);

  qp->GZwork->Norm2squared(&gap[0]); CHKERRQ(info);
  qp->TSwork->Norm2squared(&gap[1]); CHKERRQ(info);

  qp->pathnorm=sqrt( (gap[0]+gap[1]) );
  *norm=qp->pathnorm;

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_BQPIP"
static int TaoSetOptions_BQPIP(TAO_SOLVER tao, void *solver){

  TAO_BQPIP *qp = (TAO_BQPIP*)solver;
  int       info;
  TaoTruth  flg;

  TaoFunctionBegin;
  info = TaoOptionsHead("Interior point method for bound constrained quadratic optimization");CHKERRQ(info);
  info = TaoOptionInt("-predcorr","Use a predictor-corrector method","",qp->predcorr,&qp->predcorr,&flg);
  CHKERRQ(info);
  info = TaoOptionsTail();CHKERRQ(info);

  /*   info = TaoSetLinearSolverOptions(tao);CHKERRQ(info); */

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoView_BQPIP"
static int TaoView_BQPIP(TAO_SOLVER tao,void *solver){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}


/* --------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_BQPIP"
int TaoCreate_BQPIP(TAO_SOLVER tao)
{
  TAO_BQPIP *qp;
  int       info;

  TaoFunctionBegin;

  info = TaoNew(TAO_BQPIP,&qp); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_BQPIP)); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_BQPIP,(void*)qp); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_BQPIP,TaoSetDown_BQPIP); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetOptions_BQPIP); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_BQPIP); CHKERRQ(info);
  info=TaoSetTaoDualVariablesRoutine(tao,TaoGetDualVariables_BQPIP); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,100); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,500); CHKERRQ(info);
  info = TaoSetTolerances(tao,1e-12,1e-12,1e-12,0); CHKERRQ(info);
  /*
  tao->defaultmonitor     = TaoDefaultMonitor_BQPIP;
  */
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

  TaoFunctionReturn(0);
}
EXTERN_C_END
