/*$Id: s.bnls.c 1.128 02/08/16 18:01:18-05:00 benson@rockies.mcs.anl.gov $*/

#include "bnls.h"       /*I "tao_solver.h" I*/
#include "petscksp.h"
#include "petscpc.h"
#include "src/petsctao/linearsolver/taolinearsolver_petsc.h"
#include "src/petsctao/vector/taovec_petsc.h"

#include "private/kspimpl.h"
#include "private/pcimpl.h"

// Routine for BFGS preconditioner

#undef __FUNCT__
#define __FUNCT__ "bfgs_apply"
static PetscErrorCode bfgs_apply(PC pc, Vec xin, Vec xout)
{
  TaoLMVMMat *M ;
  TaoVecPetsc Xin(xin);
  TaoVecPetsc Xout(xout);
  TaoTruth info2;
  int info;

  PetscFunctionBegin;

  PetscTruth VerbosePrint = PETSC_FALSE; 
  PetscOptionsGetTruth(PETSC_NULL,"-verboseapp",&VerbosePrint,PETSC_NULL);

  info = PCShellGetContext(pc,(void**)&M); CHKERRQ(info);

  PetscScalar solnNorm,solnDot;
  info = VecNorm(xin,NORM_2,&solnNorm); CHKERRQ(info)
  info=PetscPrintf(PETSC_COMM_WORLD,"bfgs_apply: ||Xin||_2 = %22.15e\n",solnNorm);
  if(VerbosePrint) VecView(xin,0);

  info = M->Solve(&Xin, &Xout, &info2); CHKERRQ(info);

  info = VecNorm(xout,NORM_2,&solnNorm); CHKERRQ(info)
  info = VecDot(xin,xout,&solnDot); CHKERRQ(info)
  info=PetscPrintf(PETSC_COMM_WORLD,"bfgs_apply: ||Xout||_2 = %22.15e, Xin^T Xout= %22.15e\n",solnNorm,solnDot);
  if(VerbosePrint) VecView(xout,0);

  PetscFunctionReturn(0);
}

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
  
  TaoLinearSolver *tls;
  // Modify the linear solver to a conjugate gradient method
  info = TaoGetLinearSolver(tao, &tls); CHKERRQ(info);
  TaoLinearSolverPetsc *pls;
  pls  = dynamic_cast <TaoLinearSolverPetsc *> (tls);
  // set trust radius to zero 
  // PETSc ignores this case and should return the negative curvature direction
  // at its current default length
  pls->SetTrustRadius(0.0);

  if(!bnls->M) bnls->M = new TaoLMVMMat(X);
  TaoLMVMMat *M = bnls->M;
  KSP pksp = pls->GetKSP();
  // we will want to provide an initial guess in case neg curvature on the first iteration
  info = KSPSetInitialGuessNonzero(pksp,PETSC_TRUE); CHKERRQ(info);
  PC ppc;
  // Modify the preconditioner to use the bfgs approximation
  info = KSPGetPC(pksp, &ppc); CHKERRQ(info);
  PetscTruth  BFGSPreconditioner=PETSC_FALSE;// debug flag
  info = PetscOptionsGetTruth(PETSC_NULL,"-bnls_pc_bfgs",
                              &BFGSPreconditioner,PETSC_NULL); CHKERRQ(info);
  if( BFGSPreconditioner) 
    { 
     info=PetscInfo(tao,"TaoSolve_BNLS:  using bfgs preconditioner\n");
     info = KSPSetNormType(pksp, KSP_NORM_PRECONDITIONED); CHKERRQ(info);
     info = PCSetType(ppc, PCSHELL); CHKERRQ(info);
     info = PCShellSetName(ppc, "bfgs"); CHKERRQ(info);
     info = PCShellSetContext(ppc, M); CHKERRQ(info);
     info = PCShellSetApply(ppc, bfgs_apply); CHKERRQ(info);
    }
  else
    {// default to none
     info=PetscInfo(tao,"TaoSolve_BNLS:  using no preconditioner\n");
     info = PCSetType(ppc, PCNONE); CHKERRQ(info);
    }

  info = TaoComputeMeritFunctionGradient(tao,X,&f,G);CHKERRQ(info);
  info = PG->BoundGradientProjection(G,XL,X,XU);CHKERRQ(info);
  info = PG->Norm2(&gnorm); CHKERRQ(info);
  
  // Set initial scaling for the function
  if (f != 0.0) {
    info = M->SetDelta(2.0 * TaoAbsDouble(f) / (gnorm*gnorm)); CHKERRQ(info);
  }
  else {
    info = M->SetDelta(2.0 / (gnorm*gnorm)); CHKERRQ(info);
  }
  
  while (reason==TAO_CONTINUE_ITERATING){
    
    /* Project the gradient and calculate the norm */
    info = PG->BoundGradientProjection(G,XL,X,XU);CHKERRQ(info);
    info = PG->Norm2(&gnorm); CHKERRQ(info);
    
    info = M->Update(X, PG); CHKERRQ(info);

    PetscScalar ewAtol  = PetscMin(0.5,gnorm)*gnorm;
    info = KSPSetTolerances(pksp,PETSC_DEFAULT,ewAtol,
                            PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(info);
    info=PetscInfo1(tao,"TaoSolve_BNLS: gnorm =%g\n",gnorm);
    pksp->printreason = PETSC_TRUE;
    info = KSPView(pksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(info);
    M->View();

    info = TaoMonitor(tao,iter++,f,gnorm,0.0,stepsize,&reason);
    CHKERRQ(info);
    if (reason!=TAO_CONTINUE_ITERATING) break;

    info = FreeVariables->WhichEqual(PG,G); CHKERRQ(info);

    info = TaoComputeHessian(tao,X,H);CHKERRQ(info);
    
    /* Create a reduced linear system */

    info = R->SetReducedVec(G,FreeVariables);CHKERRQ(info);
    info = R->Negate();CHKERRQ(info);

    /* Use gradient as initial guess */
    PetscTruth  UseGradientIG=PETSC_FALSE;// debug flag
    info = PetscOptionsGetTruth(PETSC_NULL,"-bnls_use_gradient_ig",
                                &UseGradientIG,PETSC_NULL); CHKERRQ(info);
    if(UseGradientIG)
      info = DX->CopyFrom(G);
    else
     {
      info=PetscInfo(tao,"TaoSolve_BNLS: use bfgs init guess \n");
      info = M->Solve(G, DX, &success);
     }
    CHKERRQ(info);
    info = DXFree->SetReducedVec(DX,FreeVariables);CHKERRQ(info);
    info = DXFree->Negate(); CHKERRQ(info);
    
    info = Hsub->SetReducedMatrix(H,FreeVariables,FreeVariables);CHKERRQ(info);

    bnls->gamma_factor /= 2;
    success = TAO_FALSE;

    while (success==TAO_FALSE) {
      
      /* Approximately solve the reduced linear system */
      info = TaoPreLinearSolve(tao,Hsub);CHKERRQ(info);
      info = TaoLinearSolve(tao,Hsub,R,DXFree,&success);CHKERRQ(info);

      info = DX->SetToZero(); CHKERRQ(info);
      info = DX->ReducedXPY(DXFree,FreeVariables);CHKERRQ(info);
      info = DX->BoundGradientProjection(DX,XL,X,XU); CHKERRQ(info);
      info = DX->Dot(G,&gdx); CHKERRQ(info);

      if (gdx>=0 || success==TAO_FALSE) { /* use bfgs direction */
        info=PetscInfo1(tao,"Newton Solve Fail use BFGS direction, gdx %22.12e \n",gdx);
        info = M->Solve(G, DX, &success); CHKERRQ(info);
        info = DX->Negate(); CHKERRQ(info);
        info = DX->BoundGradientProjection(DX,XL,X,XU); CHKERRQ(info);
        // Check for success (descent direction)
        info = DX->Dot(G,&gdx); CHKERRQ(info);
        if (gdx >= 0) {
          // Step is not descent or solve was not successful
          // Use steepest descent direction (scaled)
          info=PetscInfo1(tao,"LMVM Solve Fail use steepest descent, gdx %22.12e \n",gdx);
          if (f != 0.0) {
            info = M->SetDelta(2.0 * TaoAbsDouble(f) / (gnorm*gnorm)); CHKERRQ(info);
          }
          else {
            info = M->SetDelta(2.0 / (gnorm*gnorm)); CHKERRQ(info);
          }
          info = M->Reset(); CHKERRQ(info);
          info = M->Update(X, G); CHKERRQ(info);
          info = DX->CopyFrom(G);
          info = DX->Negate(); CHKERRQ(info);
        } 
	success = TAO_TRUE;
//        bnls->gamma_factor *= 2; 
//        bnls->gamma = bnls->gamma_factor*(gnorm); 
//#if !defined(PETSC_USE_COMPLEX)
//        info=PetscInfo2(tao,"TaoSolve_NLS:  modify diagonal (assume same nonzero structure), gamma_factor=%g, gamma=%g\n",bnls->gamma_factor,bnls->gamma);
//	CHKERRQ(info);
//#else
//        info=PetscInfo3(tao,"TaoSolve_NLS:  modify diagonal (asuume same nonzero structure), gamma_factor=%g, gamma=%g, gdx %22.12e \n",
//	     bnls->gamma_factor,PetscReal(bnls->gamma),gdx);CHKERRQ(info);
//#endif
//        info = Hsub->ShiftDiagonal(bnls->gamma);CHKERRQ(info);
//        if (f != 0.0) {
//          info = M->SetDelta(2.0 * TaoAbsDouble(f) / (gnorm*gnorm)); CHKERRQ(info);
//        }
//        else {
//          info = M->SetDelta(2.0 / (gnorm*gnorm)); CHKERRQ(info);
//        }
//        info = M->Reset(); CHKERRQ(info);
//        info = M->Update(X, G); CHKERRQ(info);
//        success = TAO_FALSE;
      } else {
        info=PetscInfo1(tao,"Newton Solve is descent direction, gdx %22.12e \n",gdx);
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
  info = TaoMatDestroy(bnls->M); CHKERRQ(info);
  
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
  bnls->M=0;

  info = TaoCreateMoreThuenteBoundLineSearch(tao,0,0.9); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END
