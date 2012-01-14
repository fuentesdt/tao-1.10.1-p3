
#include "icp.h"
static int TaoSetDown_ICP(TAO_SOLVER, void*);

#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_ICP"
static int TaoSetUp_ICP(TAO_SOLVER tao, void *solver)
{
  TAO_ICP *qp =(TAO_ICP*)solver;
  int       info;
  TaoInt    n;

  TaoFunctionBegin;

  /* Set pointers to Data */
  info = TaoGetJacobian(tao,&qp->J);CHKERRQ(info);
  info = TaoGetSolution(tao,&qp->x);CHKERRQ(info);
  info = qp->x->GetDimension(&qp->n); CHKERRQ(info);

  /* Allocate some arrays */
  info = qp->x->Clone(&qp->f); CHKERRQ(info);
  info = qp->x->Clone(&qp->g2); CHKERRQ(info);
  info = qp->x->Clone(&qp->work); CHKERRQ(info);
  info = qp->x->Clone(&qp->dx); CHKERRQ(info);
  info = qp->x->Clone(&qp->JDiag); CHKERRQ(info);
  info = qp->x->Clone(&qp->diagaxpy); CHKERRQ(info);
  info = qp->x->Clone(&qp->rhs); CHKERRQ(info);
  info = qp->x->Clone(&qp->rhs2); CHKERRQ(info);
  info = qp->x->Clone(&qp->r12); CHKERRQ(info);

  info = qp->x->Clone(&qp->xl); CHKERRQ(info);
  info = qp->xl->Clone(&qp->g); CHKERRQ(info);
  info = qp->xl->Clone(&qp->dg); CHKERRQ(info);
  info = qp->xl->Clone(&qp->z); CHKERRQ(info);
  info = qp->xl->Clone(&qp->dz); CHKERRQ(info);
  info = qp->xl->Clone(&qp->gzwork); CHKERRQ(info);
  info = qp->xl->Clone(&qp->r3); CHKERRQ(info);

  info = qp->x->Clone(&qp->xu); CHKERRQ(info);
  info = qp->xu->Clone(&qp->t); CHKERRQ(info);
  info = qp->xu->Clone(&qp->dt); CHKERRQ(info);
  info = qp->xu->Clone(&qp->s); CHKERRQ(info);
  info = qp->xu->Clone(&qp->ds); CHKERRQ(info);
  info = qp->xu->Clone(&qp->tswork); CHKERRQ(info);
  info = qp->xu->Clone(&qp->r5); CHKERRQ(info);

  info = TaoSetVariableBounds(tao,qp->xl,qp->xu);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,qp->dx);CHKERRQ(info);
  info = TaoSetLagrangianGradientVector(tao,qp->g2);CHKERRQ(info);

  /* Register the Events */
  info = qp->x->GetDimension(&qp->n); CHKERRQ(info);

  qp->m=0;
  info=qp->g->GetDimension(&n); CHKERRQ(info); qp->m+=n;
  info=qp->t->GetDimension(&n); CHKERRQ(info); qp->m+=n;

  info = TaoCreateLinearSolver(tao,qp->J,300,0); CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPIPSetInitialPoint"
static int  QPIPSetInitialPoint(TAO_SOLVER tao, TAO_ICP *qp)
{
  int       info;
  double    p01=1.0e-2;
  double    mu;
  double    dd1,dd2;
  TaoVec    *g = qp->g, *z=qp->z, *s=qp->s, *t=qp->t;
  TaoVec    *gzwork=qp->gzwork, *tswork=qp->tswork;
  TaoVec    *rhs=qp->rhs, *r12=qp->r12;
  TaoVec    *x = qp->x, *xl=qp->xl, *xu=qp->xu, *f=qp->f;
  TaoVec    *work=qp->work;

  TaoFunctionBegin;

  /* Initialize Primal Vectors */
  info = t->Waxpby(1.0,xu,-1.0,x);CHKERRQ(info);
  info = g->Waxpby(1.0,x,-1.0,xl);CHKERRQ(info);

  info = gzwork->SetToConstant(p01);CHKERRQ(info);
  info = tswork->SetToConstant(p01);CHKERRQ(info);

  info = g->PointwiseMaximum(g,gzwork); CHKERRQ(info);
  info = t->PointwiseMaximum(t,tswork); CHKERRQ(info);

  /* Initialize Dual Variable Vectors */

  info = z->CopyFrom(g); CHKERRQ(info);
  info = z->Reciprocal(); CHKERRQ(info);

  info = s->CopyFrom(t); CHKERRQ(info);
  info = s->Reciprocal(); CHKERRQ(info);
  
  info = r12->Waxpby(1.0,s,-1.0,z);CHKERRQ(info);
  info = rhs->CopyFrom(f);CHKERRQ(info);
  
  info = rhs->AbsoluteValue(); CHKERRQ(info);
  info = r12->AbsoluteValue(); CHKERRQ(info);
  info = work->SetToConstant(p01);CHKERRQ(info);
  info = rhs->PointwiseMaximum(rhs,work); CHKERRQ(info);
  info = r12->PointwiseMaximum(r12,work); CHKERRQ(info);
  
  info = rhs->PointwiseDivide(rhs,r12); CHKERRQ(info);
  
  info = r12->Norm1(&qp->dinfeas); CHKERRQ(info);
  mu = (qp->dinfeas)/(qp->n);
  
  info = s->Scale(mu); CHKERRQ(info);
  info = z->Scale(mu); CHKERRQ(info);
  
  info = QPIPComputeResidual(qp); CHKERRQ(info);   
  info = TAOComputeNormFromCentralPath_ICP(tao,&dd2); CHKERRQ(info);
  dd1=TaoMax(qp->pathnorm,1.0e-2);
  while ( (qp->dinfeas+qp->pinfeas)/(qp->m+qp->n) >= qp->mu ){
    
    info=g->AddConstant(dd1); CHKERRQ(info);
    info=z->AddConstant(dd1); CHKERRQ(info);
    info=s->AddConstant(dd1); CHKERRQ(info);
    info=t->AddConstant(dd1); CHKERRQ(info);

    info = QPIPComputeResidual(qp); CHKERRQ(info);   
    info = TAOComputeNormFromCentralPath_ICP(tao,&dd2); CHKERRQ(info);
    
    qp->dobj = qp->pobj - qp->gap;
    qp->rgap=qp->gap/( TaoAbsScalar(qp->dobj) + TaoAbsScalar(qp->pobj) + 1.0 );
    dd1=dd1*2;
  }

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_ICP"
static int TaoSetDown_ICP(TAO_SOLVER tao, void*solver)
{
  TAO_ICP *qp = (TAO_ICP*)solver;
  int info;

  /* Free allocated memory in GPCG structure */
  TaoFunctionBegin;

  info=TaoVecDestroy(qp->g2);CHKERRQ(info);
  info=TaoVecDestroy(qp->g);CHKERRQ(info);
  info=TaoVecDestroy(qp->dg);CHKERRQ(info);
  info=TaoVecDestroy(qp->z);CHKERRQ(info);
  info=TaoVecDestroy(qp->dz);CHKERRQ(info);
  info=TaoVecDestroy(qp->gzwork);CHKERRQ(info);
  info=TaoVecDestroy(qp->r3);CHKERRQ(info);
  info=TaoVecDestroy(qp->xl);CHKERRQ(info);
  
  info=TaoVecDestroy(qp->s);CHKERRQ(info);
  info=TaoVecDestroy(qp->ds);CHKERRQ(info);
  info=TaoVecDestroy(qp->t);CHKERRQ(info);
  info=TaoVecDestroy(qp->dt);CHKERRQ(info);
  info=TaoVecDestroy(qp->tswork);CHKERRQ(info);
  info=TaoVecDestroy(qp->r5);CHKERRQ(info);
  info=TaoVecDestroy(qp->xu);CHKERRQ(info);
  
  info=TaoVecDestroy(qp->f);CHKERRQ(info);
  info=TaoVecDestroy(qp->JDiag);CHKERRQ(info);
  info=TaoVecDestroy(qp->work);CHKERRQ(info);
  info=TaoVecDestroy(qp->diagaxpy);CHKERRQ(info);
  info=TaoVecDestroy(qp->rhs);CHKERRQ(info);
  info=TaoVecDestroy(qp->rhs2);CHKERRQ(info);
  info=TaoVecDestroy(qp->dx);CHKERRQ(info);
  info=TaoVecDestroy(qp->r12);CHKERRQ(info);
  
  info = TaoDestroyLinearSolver(tao);CHKERRQ(info);

  info = TaoSetVariableBounds(tao,0,0);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,0);CHKERRQ(info);
  info = TaoSetLagrangianGradientVector(tao,qp->g2);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_ICP"
static int TaoSolve_ICP(TAO_SOLVER tao, void *solver)
{
  TAO_ICP   *qp = (TAO_ICP*)solver;
  TaoVec    *x = qp->x, *xl=qp->xl, *xu=qp->xu, *f=qp->f;
  TaoVec    *g = qp->g, *z=qp->z, *s=qp->s, *t=qp->t;
  TaoVec    *dg=qp->dg, *dz=qp->dz, *ds=qp->ds, *dt=qp->dt;
  TaoVec    *gzwork=qp->gzwork, *tswork=qp->tswork;
  TaoVec    *rhs=qp->rhs, *rhs2=qp->rhs2;
  TaoVec    *diagaxpy=qp->diagaxpy, *dx=qp->dx;
  TaoVec    *r3=qp->r3,*r5=qp->r5;
  TaoVec    *JDiag=qp->JDiag;
  TaoMat    *J;
  int       info;
  TaoInt    iter=0;
  TaoInt    n;
  double    d1,ksptol;
  double    sigmamu;// sigma;
  double    dstep,pstep,step=0;
  double    gap[4];
  TaoTruth  kspsuccess;
  TaoTerminateReason reason;
  TaoLinearSolver *qpksp;

  TaoFunctionBegin;

  info = TaoGetSolution(tao,&x);CHKERRQ(info); qp->x=x;
  info = TaoGetJacobian(tao,&J);CHKERRQ(info); qp->J=J;

  info = x->Median(xl,x,xu); CHKERRQ(info);
  info = TaoEvaluateVariableBounds(tao,xl,xu); CHKERRQ(info);
  info = TaoComputeConstraints(tao, x, f); CHKERRQ(info);
 
  info = x->GetDimension(&n);CHKERRQ(info);

  info = QPIPSetInitialPoint(tao,qp); CHKERRQ(info);
  
  /* Enter main loop */
  while (1){

    info = QPIPComputeResidual(qp); CHKERRQ(info);
    info = TAOComputeNormFromCentralPath_ICP(tao,&d1); CHKERRQ(info);
    
    /* Check Stopping Condition      */
    info=TaoMonitor(tao,iter++,qp->gap,(qp->gap+qp->dinfeas),qp->pinfeas,
		    step,&reason); CHKERRQ(info);
    if (reason != TAO_CONTINUE_ITERATING) break;

    info = rhs->ScaleCopyFrom(-1.0,f); CHKERRQ(info);
    
    info = diagaxpy->SetToZero(); CHKERRQ(info);
    
    info = gzwork->PointwiseDivide(z,g); CHKERRQ(info);
    info = diagaxpy->Axpy(1.0,gzwork); CHKERRQ(info);
    info = gzwork->PointwiseMultiply(gzwork,r3); CHKERRQ(info);
    info = rhs->Axpy(-1.0,gzwork); CHKERRQ(info);
    
    info = tswork->PointwiseDivide(s,t); CHKERRQ(info);
    info = diagaxpy->Axpy(1.0,tswork); CHKERRQ(info);
    info = tswork->PointwiseMultiply(tswork,r5); CHKERRQ(info);
    info = rhs->Axpy(-1.0,tswork); CHKERRQ(info);
    /*    
    info = dz->SetToZero(); CHKERRQ(info);
    info = dg->SetToZero(); CHKERRQ(info);
    info = ds->SetToZero(); CHKERRQ(info);
    info = dt->SetToZero(); CHKERRQ(info);
    */
    /* Add centrality part */
    if (iter > 0 && (qp->rnorm>5*qp->mu || d1*d1>qp->m*qp->mu*qp->mu) ) {
      sigmamu=qp->mu;//sigma=1.0;
      sigmamu=0;//sigma=0.0;
    } else {
      sigmamu=0;//sigma=0.0;
    }
    sigmamu=0.1*qp->mu;//sigma=0.1;
    sigmamu=0.0*qp->mu;//sigma=0.0;
    info = dz->SetToConstant(sigmamu); CHKERRQ(info);
    info = ds->SetToConstant(sigmamu); CHKERRQ(info);
    //    printf("TargetMu: %4.2e\n",sigmamu);
    if (sigmamu !=0){
      info = dz->PointwiseDivide(dz,g); CHKERRQ(info);
      info = ds->PointwiseDivide(ds,t); CHKERRQ(info);
      info = rhs->Axpy(1.0,dz); CHKERRQ(info);
      info = rhs->Axpy(-1.0,ds); CHKERRQ(info);
    } else {
      //      info = qp->RHS2->SetToZero(); CHKERRQ(info);
    }

    qp->usedcorrector=TAO_FALSE;
    qp->sigmamu=sigmamu;



    /*  Determine the solving tolerance */
    ksptol = qp->mu/10.0;
    ksptol = TaoMin(ksptol,0.00000000001);
    info = TaoGetLinearSolver(tao,&qpksp); CHKERRQ(info);
    /*
    info = qpksp->SetTolerances(ksptol,1e-30,1e30,2*n); CHKERRQ(info);
    */

    info = TaoComputeJacobian(tao,x,J); CHKERRQ(info);
    info = J->GetDiagonal(JDiag); CHKERRQ(info);
    info = J->AddDiagonal(diagaxpy); CHKERRQ(info);

    info = TaoPreLinearSolve(tao,J);CHKERRQ(info);
    info = TaoLinearSolve(tao,J,rhs,dx,&kspsuccess);CHKERRQ(info);

    info = QPComputeStepDirection(qp); CHKERRQ(info);
    info = QPStepLength(qp);  CHKERRQ(info);

    /* Calculate New Residual R1 in Work vector */
    info=dg->Dot(dz, gap); CHKERRQ(info);
    info=dt->Dot(ds, gap+1); CHKERRQ(info);
 
    qp->rnorm=(qp->dinfeas+qp->psteplength*qp->pinfeas)/(qp->m+qp->n);
    pstep = qp->psteplength; dstep = qp->dsteplength;    
    step = TaoMin(qp->psteplength,qp->dsteplength);
    //    printf("Step: %4.2f\n",step);
    step=1-step;
    sigmamu = step*step*step*qp->mu;
    //    sigmamu = 0;
    qp->sigmamu=sigmamu;
    /*
    ( pstep*pstep*(gap[0]+gap[1]) +
	       (1 - pstep + pstep*sigma)*qp->gap  )/qp->m;
    */
    if (1==1 && /* && qp->predcorr &&*/ step < 1.99){
      
      /*
      if (sigmamu < qp->mu){ 
	sigmamu=sigmamu/qp->mu;
	sigmamu=sigmamu*sigmamu*sigmamu;
      } else {sigmamu = 1.0;}
      sigmamu = 0.0*qp->mu;
      */
      qp->usedcorrector=TAO_TRUE;
      /* Compute Corrector Step */
      info = gzwork->PointwiseMultiply(dg,dz); CHKERRQ(info);
      info = gzwork->Negate(); CHKERRQ(info);
      info = gzwork->AddConstant(sigmamu); CHKERRQ(info);
      info = gzwork->PointwiseDivide(gzwork,g); CHKERRQ(info);

      info = tswork->PointwiseMultiply(ds,dt); CHKERRQ(info);
      info = tswork->Negate(); CHKERRQ(info);
      info = tswork->AddConstant(sigmamu); CHKERRQ(info);
      info = tswork->PointwiseDivide(tswork,t); CHKERRQ(info);
      
      info = rhs2->Waxpby(1.0,gzwork,-1.0,tswork); CHKERRQ(info);
      info = rhs2->Axpy(1.0,rhs); CHKERRQ(info);

      info = TaoLinearSolve(tao,J,rhs2,dx,&kspsuccess);CHKERRQ(info);
      
      info = QPComputeStepDirection(qp); CHKERRQ(info);
      info = QPStepLength(qp); CHKERRQ(info);

    }  /* End Corrector step */


    /* Take the step */
    pstep = qp->psteplength; dstep = qp->dsteplength;

    info = z->Axpy(dstep,dz); CHKERRQ(info);
    info = s->Axpy(dstep,ds); CHKERRQ(info);
    info = x->Axpy(dstep,dx); CHKERRQ(info);
    info = g->Axpy(pstep,dg); CHKERRQ(info);
    info = t->Axpy(pstep,dt); CHKERRQ(info);
    
    info = TaoComputeConstraints(tao, x, f); CHKERRQ(info);

  }  /* END MAIN LOOP  */

  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPComputeStepDirection"
static int QPComputeStepDirection(TAO_ICP *qp)
{
  int info;
  double gap1,gap2;
  TaoVec    *g = qp->g, *z=qp->z, *s=qp->s, *t=qp->t;
  TaoVec    *dg=qp->dg, *dz=qp->dz, *ds=qp->ds, *dt=qp->dt;
  TaoVec    *gzwork=qp->gzwork, *tswork=qp->tswork;
  TaoVec    *dx=qp->dx;
  TaoVec    *r3=qp->r3,*r5=qp->r5;

  TaoFunctionBegin;

  if (qp->usedcorrector==TAO_TRUE){
    info = dz->PointwiseMultiply(dz,dg); CHKERRQ(info);
    info = dz->PointwiseDivide(dz,g); CHKERRQ(info);
    info = dz->Negate(); CHKERRQ(info);
    
    info = ds->PointwiseMultiply(ds,dt); CHKERRQ(info);
    info = ds->PointwiseDivide(ds,t); CHKERRQ(info);
    info = ds->Negate(); CHKERRQ(info);
  } 
  /*
    info = dz->SetToZero(); CHKERRQ(info);
    info = ds->SetToZero(); CHKERRQ(info);
  */
   /* Calculate DG */
  info = dg->Waxpby(1.0,dx,1.0,r3);CHKERRQ(info);

  /* Calculate DT */
  info = dt->Waxpby(-1.0,dx,-1.0,r5);CHKERRQ(info);

  /* Calculate DZ */
  info = gzwork->PointwiseMultiply(dg,z); CHKERRQ(info);
  info = gzwork->Negate(); CHKERRQ(info);
  info = gzwork->AddConstant(qp->sigmamu); CHKERRQ(info);
  info = gzwork->PointwiseDivide(gzwork,g); CHKERRQ(info);
  info = gzwork->Axpy(-1.0,z);CHKERRQ(info);

  info = dz->Axpy(1.0,gzwork);CHKERRQ(info);
  //  info = qp->DZ->CopyFrom(qp->GZwork);CHKERRQ(info);

  /* Calculate DS */
  info = tswork->PointwiseMultiply(dt,s); CHKERRQ(info);
  info = tswork->Negate(); CHKERRQ(info);
  info = tswork->AddConstant(qp->sigmamu); CHKERRQ(info);
  info = tswork->PointwiseDivide(tswork,t); CHKERRQ(info);
  info = tswork->Axpy(-1.0,s);CHKERRQ(info);
  info = ds->Axpy(1.0,tswork);CHKERRQ(info);
  //  info = qp->DS->CopyFrom(qp->TSwork);CHKERRQ(info);

  info = ds->Dot(dt,&gap1); CHKERRQ(info);
  info = dg->Dot(dz,&gap2); CHKERRQ(info);
  //  printf("Gap Reduce: %4.2e \n",gap1+gap2);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPIPComputeResidual"
static int QPIPComputeResidual(TAO_ICP *qp)
{
  int info;
  double gap1,gap2,dtmp = 1.0 - qp->psteplength;
  TaoVec    *x = qp->x, *xl=qp->xl, *xu=qp->xu, *f=qp->f;
  TaoVec    *g = qp->g, *z=qp->z, *s=qp->s, *t=qp->t;
  TaoVec    *rhs=qp->rhs;
  TaoVec    *r12=qp->r12,*r3=qp->r3,*r5=qp->r5;

  TaoFunctionBegin;

  /* Compute R3 and R5 */

  if (0==1){

    info = r3->Scale(dtmp);
    info = r5->Scale(dtmp);
    qp->pinfeas=dtmp*qp->pinfeas;

  } else {

    info = r3->Waxpby(1.0,x,-1.0,xl);CHKERRQ(info);
    info = r3->Axpy(-1.0,g);CHKERRQ(info);

    info = r5->Waxpby(1.0,x,-1.0,xu);CHKERRQ(info);
    info = r5->Axpy(1.0,t);CHKERRQ(info);

    info = r3->NormInfinity(&gap1);CHKERRQ(info);
    info = r5->NormInfinity(&gap2);CHKERRQ(info);

    qp->pinfeas=TaoMax(gap1,gap2);

  }

  info = r12->Waxpby(1.0,s,-1.0,z);CHKERRQ(info);

  info = rhs->CopyFrom(f);CHKERRQ(info);
  info = rhs->Negate();CHKERRQ(info);

  info = r12->Axpy(-1.0,rhs);CHKERRQ(info);

  info = r12->NormInfinity(&qp->dinfeas); CHKERRQ(info);
  qp->rnorm=(qp->dinfeas+qp->pinfeas)/(qp->m+qp->n);
  qp->rnorm=TaoMax(qp->dinfeas,qp->pinfeas);
  //  printf("Dual Infeas: %4.2e,  Primal Infeas: %4.2e\n",qp->dinfeas,qp->pinfeas);
  
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "QPStepLength"
static int QPStepLength(TAO_ICP *qp)
{
  double tstep1,tstep2,tstep3,tstep4,tstep;
  int info;
  TaoVec    *g = qp->g, *z=qp->z, *s=qp->s, *t=qp->t;
  TaoVec    *dg=qp->dg, *dz=qp->dz, *ds=qp->ds, *dt=qp->dt;

  TaoFunctionBegin;
  /* Compute stepsize to the boundary */

  info = g->StepMax(dg,&tstep1); CHKERRQ(info);
  info = t->StepMax(dt,&tstep2); CHKERRQ(info);
  info = s->StepMax(ds,&tstep3); CHKERRQ(info);
  info = z->StepMax(dz,&tstep4); CHKERRQ(info);

  tstep = TaoMin(tstep1,tstep2);
  qp->psteplength = TaoMin(0.999*tstep,1.0);

  tstep = TaoMin(tstep3,tstep4);
  qp->dsteplength = TaoMin(0.999*tstep,1.0);

  qp->psteplength = TaoMin(qp->psteplength,qp->dsteplength);
  qp->dsteplength = qp->psteplength;
  //  printf("StepSize: %4.2e\n",qp->dsteplength);

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TAOComputeNormFromCentralPath_ICP"
int TAOComputeNormFromCentralPath_ICP(TAO_SOLVER tao, double *norm)
{
  TAO_ICP *qp;
  int       info;
  double    gap[2],mu[2], nmu;
  
  TaoFunctionBegin;
  info = TaoGetSolverContext(tao,"tao_icp",(void**)&qp); CHKERRQ(info);
  TaoVec    *g = qp->g, *z=qp->z, *s=qp->s, *t=qp->t;
  TaoVec    *gzwork=qp->gzwork, *tswork=qp->tswork;

  info = gzwork->PointwiseMultiply(g,z); CHKERRQ(info);
  info = tswork->PointwiseMultiply(t,s); CHKERRQ(info);
  info = tswork->Norm1(&mu[0]); CHKERRQ(info);
  info = gzwork->Norm1(&mu[1]); CHKERRQ(info);

  nmu=-(mu[0]+mu[1])/qp->m;
  qp->gap=mu[0]+mu[1];
  qp->mu=-nmu;

  info = gzwork->AddConstant(nmu); CHKERRQ(info);
  info = tswork->AddConstant(nmu); CHKERRQ(info);

  info = gzwork->NormInfinity(&gap[0]); CHKERRQ(info);
  info = tswork->NormInfinity(&gap[1]); CHKERRQ(info);

  qp->pathnorm=TaoMax(gap[0],gap[1]);
  *norm=qp->pathnorm;
  //  printf("Gap: %4.2e,  Mu: %4.2e,  Norm: %4.2e\n",qp->gap,-nmu,qp->pathnorm);

  //  printf("Gap: %4.2e,  Mu: %4.2e,  Norm: %4.2e\n",qp->gap,qp->mu,qp->pathnorm);
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_ICP"
static int TaoSetOptions_ICP(TAO_SOLVER tao, void *solver){

  TAO_ICP *qp = (TAO_ICP*)solver;
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
#define __FUNCT__ "TaoView_ICP"
static int TaoView_ICP(TAO_SOLVER tao,void *solver){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}


/* --------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_ICP"
int TaoCreate_ICP(TAO_SOLVER tao)
{
  TAO_ICP *qp;
  int       info;

  TaoFunctionBegin;

  info = TaoNew(TAO_ICP,&qp); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_ICP)); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_ICP,(void*)qp); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_ICP,TaoSetDown_ICP); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetOptions_ICP); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_ICP); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,100); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,500); CHKERRQ(info);
  info = TaoSetTolerances(tao,1e-10,1e-10,1e-12,0); CHKERRQ(info);
  /*
  tao->defaultmonitor     = TaoDefaultMonitor_ICP;
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
