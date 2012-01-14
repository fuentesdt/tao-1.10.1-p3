#include "src/complementarity/impls/ssls/ssls.h"
// #include "src/tao_impl.h"

int Tao_RSCS_FunctionGradient(TAO_SOLVER, TaoVec*, double*, TaoVec*, void *);
int Tao_RSCS_Function(TAO_SOLVER, TaoVec *, double *, void *);
static int TaoSolve_RSCS(TAO_SOLVER tao, void*solver);

typedef struct {
  TaoVec *f;
  TaoMat *J;
} FGMeritCtx;


#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_RSCS"
static int TaoSolve_RSCS(TAO_SOLVER tao, void*solver){

  TAO_SSLS *asls = (TAO_SSLS *)solver;
  FGMeritCtx  meritctx;
  TaoTerminateReason reason=TAO_CONTINUE_ITERATING;
  TaoVec *x, *l, *u, *ff, *d, *w, *g;
  TaoMat *J,*Jsub;
  TaoVec *dxfree,*r1;
  TaoIndexSet *FreeVariableSet;
  int  info;
  TaoInt lsflag,iter=0;
  double ndpsi,stepsize=1.0;
  //  double psi,psi1,psi2;
  TaoTruth success;
  double fff,fff_full,gdx;
  double gamma = 0.0, gamma_factor =1.0e-12;
  //  double fref;

  TaoFunctionBegin;

  ff = asls->dpsi;
  d = asls->d;
  g = asls->t1;
  w = asls->w;

  /*
  f = asls->f;
  ff = asls->ff;
  t2 = asls->t2;
  da = asls->da;
  db = asls->db;
  */
  /* Check if upper bound greater than lower bound. */
  info = TaoGetSolution(tao, &x); CHKERRQ(info);
  info = TaoGetVariableBounds(tao, &l, &u); CHKERRQ(info);
  info = TaoEvaluateVariableBounds(tao,l,u); CHKERRQ(info);
  info = TaoGetJacobian(tao, &J); CHKERRQ(info);
  meritctx.J=J;
  meritctx.f=ff;
  info = TaoSetMeritFunction(tao, Tao_RSCS_Function, Tao_RSCS_FunctionGradient,
                             TAO_NULL, TAO_NULL, TAO_NULL, (void*)&meritctx ); CHKERRQ(info);

  /*   Project the current point onto the feasible set */
  info = x->Median(l,x,u); CHKERRQ(info);
  
  info = TaoComputeMeritFunctionGradient(tao, x, &fff, g); CHKERRQ(info);

  info = x->CreateIndexSet(&FreeVariableSet); CHKERRQ(info);
  info = J->CreateReducedMatrix(FreeVariableSet,FreeVariableSet,&Jsub); CHKERRQ(info);
  info = x->Clone(&dxfree); CHKERRQ(info);
  info = x->Clone(&r1); CHKERRQ(info);

  while (reason==TAO_CONTINUE_ITERATING){
    
    /* Project the gradient and calculate the norm */
    info = w->BoundGradientProjection(ff,l,x,u);CHKERRQ(info);
    info = w->Norm2(&ndpsi); CHKERRQ(info);

    info = TaoMonitor(tao,iter++,fff,ndpsi,0.0,stepsize,&reason);
    CHKERRQ(info);

    if (reason!=TAO_CONTINUE_ITERATING) break;

    info = FreeVariableSet->WhichEqual(w,ff); CHKERRQ(info);

    /* Create a reduced linear system */
    info = r1->SetReducedVec(ff,FreeVariableSet);CHKERRQ(info);
    info = r1->Negate();CHKERRQ(info);

    info = dxfree->SetReducedVec(d,FreeVariableSet);CHKERRQ(info);
    info = dxfree->SetToZero(); CHKERRQ(info);
    
    info = TaoComputeJacobian(tao,x,J);CHKERRQ(info);    
    info = Jsub->SetReducedMatrix(J,FreeVariableSet,FreeVariableSet);CHKERRQ(info);
    
    success = TAO_FALSE;
    gamma = gamma_factor*(ndpsi); 
    while (success==TAO_FALSE) {
      
      /* Approximately solve the reduced linear system */
      info = TaoPreLinearSolve(tao,Jsub);CHKERRQ(info);
      info = TaoLinearSolve(tao,Jsub,r1,dxfree,&success);CHKERRQ(info);
    
      info = d->SetToZero(); CHKERRQ(info);
      info = d->ReducedXPY(dxfree,FreeVariableSet);CHKERRQ(info);
      
      info = d->Dot(ff,&gdx); CHKERRQ(info);
      

      if (success==TAO_FALSE) { /* Modify diagonal of Hessian if not a descent direction */

	info = Jsub->SetReducedMatrix(J,FreeVariableSet,FreeVariableSet);CHKERRQ(info);
        gamma *=10; 
	//	  printf("Shift diagonal: %4.2e\n",gamma);
        info = PetscInfo2(tao,"TaoSolve_NLS:  modify diagonal (asuume same nonzero structure), gamma_factor=%g, gamma=%g\n",gamma_factor,gamma);CHKERRQ(info);
        info = Jsub->ShiftDiagonal(gamma);CHKERRQ(info);
	success = TAO_FALSE;
	
      } else {
	success = TAO_TRUE;
      }
      
    }

    //    fref=fff;
    info = g->ScaleCopyFrom(-1.0,d);CHKERRQ(info);
    stepsize=1.0;	
    info = TaoLineSearchApply(tao,x,g,d,w,
			      &fff,&fff_full, &stepsize,&lsflag);
    CHKERRQ(info);


    if (lsflag!=0){
      int attempts = 0;

      printf("Try Again \n");
      info = FreeVariableSet->WhichBetween(l,x,u); CHKERRQ(info);
      
      /* Create a reduced linear system */
      info = r1->SetReducedVec(ff,FreeVariableSet);CHKERRQ(info);
      info = r1->Negate();CHKERRQ(info);
      
      info = dxfree->SetReducedVec(d,FreeVariableSet);CHKERRQ(info);
      info = dxfree->SetToZero(); CHKERRQ(info);
      
      info = TaoComputeJacobian(tao,x,J);CHKERRQ(info);    
      info = Jsub->SetReducedMatrix(J,FreeVariableSet,FreeVariableSet);CHKERRQ(info);
      
      success = TAO_FALSE;
      gamma = gamma_factor*(ndpsi); 
      while (success==TAO_FALSE && attempts < 10) {
	
	/* Approximately solve the reduced linear system */
	info = TaoPreLinearSolve(tao,Jsub);CHKERRQ(info);
	info = TaoLinearSolve(tao,Jsub,r1,dxfree,&success);CHKERRQ(info);
	
	info = d->SetToZero(); CHKERRQ(info);
	info = d->ReducedXPY(dxfree,FreeVariableSet);CHKERRQ(info);
	
	info = d->Dot(ff,&gdx); CHKERRQ(info);
	
	if (success==TAO_FALSE) { /* Modify diagonal of Hessian if not a descent direction */
	  
	  info = Jsub->SetReducedMatrix(J,FreeVariableSet,FreeVariableSet);CHKERRQ(info);
	  gamma *= 10; 
	  //	  printf("Shift diagonal: %4.2e\n",gamma);
	  info = PetscInfo2(tao,"TaoSolve_NLS:  modify diagonal (asuume same nonzero structure), gamma_factor=%g, gamma=%g\n",
			       gamma_factor,gamma); CHKERRQ(info);
	  info = Jsub->ShiftDiagonal(gamma);CHKERRQ(info);
	  success = TAO_FALSE;
	  
	} else {
	  success = TAO_TRUE;
	}
        ++attempts;
      }
      
      //      fref=fff;
      info = g->ScaleCopyFrom(-1.0,d);CHKERRQ(info);
      stepsize=1.0;	
      info = TaoLineSearchApply(tao,x,g,d,w,
				&fff,&fff_full, &stepsize,&lsflag);
      CHKERRQ(info);
      
    }

    if ( iter>=40 && (iter%10==0) ){
      printf("Steepest Descent \n");
      info = d->ScaleCopyFrom(-1.0,ff);CHKERRQ(info);
      info = g->CopyFrom(ff);CHKERRQ(info);
      stepsize=1.0;	
      info = TaoLineSearchApply(tao,x,g,d,w,
				&fff,&fff_full, &stepsize,&lsflag);
      CHKERRQ(info);
      stepsize = 1.0e-6;
    }

    if (lsflag!=0 ){
      printf("Steepest Descent \n");
      info = d->ScaleCopyFrom(-1.0,ff);CHKERRQ(info);
      info = g->CopyFrom(ff);CHKERRQ(info);
      stepsize=1.0;	
      info = TaoLineSearchApply(tao,x,g,d,w,
				&fff,&fff_full, &stepsize,&lsflag);
      CHKERRQ(info);
      stepsize = 1.0e-6;
      
    }

    printf("Stepsize: %4.2e\n",stepsize);
  }  /* END MAIN LOOP  */


  info = TaoMatDestroy(Jsub);CHKERRQ(info);
  info = TaoVecDestroy(dxfree);CHKERRQ(info);
  info = TaoVecDestroy(r1);CHKERRQ(info);
  info = TaoIndexSetDestroy(FreeVariableSet);CHKERRQ(info);

  TaoFunctionReturn(0);
}


/* ---------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_RSCS"
int TaoCreate_RSCS(TAO_SOLVER tao)
{
  TAO_SSLS *asls;
  int        info;

  TaoFunctionBegin;

  info = TaoNew(TAO_SSLS, &asls); CHKERRQ(info);
  info = PetscLogObjectMemory(tao, sizeof(TAO_SSLS)); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_RSCS,(void*)asls); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_SSLS,TaoSetDown_SSLS); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetOptions_SSLS); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_SSLS); CHKERRQ(info);

  // info = TaoCreateMoreThuenteBoundLineSearch(tao,0,0.9); CHKERRQ(info);
  info = TaoCreateNDProjectedArmijoLineSearch(tao); CHKERRQ(info);
  //   info = TaoCreateProjectedArmijoLineSearch(tao); CHKERRQ(info);
  info = TaoSetMeritFunction(tao, Tao_RSCS_Function, Tao_RSCS_FunctionGradient,
                             TAO_NULL, TAO_NULL, TAO_NULL, asls); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,2000); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,4000); CHKERRQ(info);

  info = TaoSetTolerances(tao,0,0,0,0); CHKERRQ(info);
  info = TaoSetGradientTolerances(tao,1.0e-16,0.0,0.0); CHKERRQ(info);
  info = TaoSetFunctionLowerBound(tao,1.0e-8); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "Tao_RSCS_Function"
int Tao_RSCS_Function(TAO_SOLVER tao, TaoVec *X, double *fcn, void *solver) 
{
  int info;
  double ndpsi;
  TaoVec *TT,*XL,*XU;
  FGMeritCtx*ctx = (FGMeritCtx*)solver;

  TaoFunctionBegin;
  info = TaoGetVariableBounds(tao, &XL, &XU); CHKERRQ(info);
  info = X->Clone(&TT); CHKERRQ(info);
  info = TaoComputeConstraints(tao, X, ctx->f); CHKERRQ(info);
  info = TT->Fischer(X, ctx->f, XL, XU); CHKERRQ(info);
  //  info = TT->BoundGradientProjection(ctx->f,XL,X,XU);CHKERRQ(info);
  info = TT->Norm2(&ndpsi); CHKERRQ(info);
  *fcn=ndpsi;
  info = PetscInfo1(tao,"RSCS Function value: %4.2e\n",*fcn); CHKERRQ(info);
  info = TaoVecDestroy(TT); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "Tao_RSCS_FunctionGradient"
int Tao_RSCS_FunctionGradient(TAO_SOLVER tao, TaoVec *X, double *fcn, 
                              TaoVec *G, void *solver)
{
  int info;
  double ndpsi;
  TaoVec *XL,*XU,*TT,*f;
  TaoMat *J;
  FGMeritCtx*ctx = (FGMeritCtx*)solver;

  TaoFunctionBegin;
  J=ctx->J;
  f=ctx->f;
  info = TaoGetVariableBounds(tao, &XL, &XU); CHKERRQ(info);
  info = TaoComputeConstraints(tao, X, f); CHKERRQ(info);
  //  info = TaoComputeJacobian(tao,X,J); CHKERRQ(info);

  if (0==1){
    info = G->BoundGradientProjection(f,XL,X,XU);CHKERRQ(info);
    info = G->Norm2(&ndpsi); CHKERRQ(info);
    info = G->CopyFrom(f);CHKERRQ(info);
    *fcn=ndpsi;
  } else if (2==2){
    info = G->Fischer(X, f, XL, XU); CHKERRQ(info);
    info = G->Norm2(fcn); CHKERRQ(info);
    info = G->CopyFrom(f);CHKERRQ(info);
  } else {
    info = G->Clone(&TT); CHKERRQ(info);
    info = TT->BoundGradientProjection(f,XL,X,XU);CHKERRQ(info);
    info = TT->Norm2(&ndpsi); CHKERRQ(info);
    *fcn=ndpsi*ndpsi/2;
    *fcn=ndpsi;
    info = J->MultiplyTranspose(TT,G);CHKERRQ(info);
    info = TaoVecDestroy(TT); CHKERRQ(info);
  }

  info = PetscInfo1(tao,"RSCS Function value: %4.2e\n",*fcn); CHKERRQ(info);

  TaoFunctionReturn(0);
}
