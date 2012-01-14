#include "tao.h"
#ifdef TAO_USE_PETSC
#include "src/petsctao/application/taoapp/taoapp_petsc.h"

typedef struct {
  PetscTruth check_gradient;
  PetscTruth check_hessian;
  PetscTruth complete_print;
} TAO_FD;

#undef __FUNCT__
#define __FUNCT__ "TaoSolve_FD"
static int TaoSolve_FD(TAO_SOLVER tao, void *solver)
{
  Mat A,Apre,B;
  MatStructure flg;
  Vec            x,g1,g2;
  PetscInt       i;
  PetscReal      hcnorm,diffnorm;
  int info;
  TAO_FD *fd       = (TAO_FD*)solver;
  TaoApplication *tapp;
  TaoPetscApplication *tpapp;
  TAO_APPLICATION app;
  MPI_Comm comm, gcomm;
  TaoFunctionBegin;
  
  info = TaoGetApplication(tao,&tapp); CHKERRQ(info);
  tpapp = dynamic_cast<TaoPetscApplication*>(tapp);
  if (!tpapp) {
    SETERRQ(1,"tao_fd_check solver only works with PETSc applications");
  }
  app = tpapp->papp;
  info = tpapp->GetComm(&comm); CHKERRQ(info);
  info = TaoAppGetSolutionVec(app,&x); CHKERRQ(info);
  if (fd->check_gradient) {
    info = PetscPrintf(comm,"Testing hand-coded gradient, if the ratio ||fd - hc|| / ||hc|| is\n"); CHKERRQ(info);
    info = PetscPrintf(comm,"0 (1.e-8), the hand-coded gradient is probably correct.\n"); CHKERRQ(info);

    if (!fd->complete_print) {
      info = PetscPrintf(comm,"Run with -tao_test_display to show difference\n");CHKERRQ(info);
      info = PetscPrintf(comm,"between hand-coded and finite difference gradient.\n");CHKERRQ(info);
    }
    
    info = VecDuplicate(x,&g1); CHKERRQ(info);
    info = VecDuplicate(x,&g2); CHKERRQ(info);
    for (i=0; i<3; i++) {  
      if (i == 1) {info = VecSet(x,-1.0);CHKERRQ(info);}
      else if (i == 2) {info = VecSet(x,1.0);CHKERRQ(info);}
 
      /* compute both versions of gradient */
      info = TaoAppComputeGradient(app,x,g1); CHKERRQ(info);
      info = TaoAppDefaultComputeGradient(app,x,g2,PETSC_NULL); CHKERRQ(info);
      if (fd->complete_print) {
	info = PetscPrintf(comm,"Finite difference gradient\n"); CHKERRQ(info);
	info = PetscObjectGetComm((PetscObject)g2,&gcomm); CHKERRQ(info);
	info = VecView(g2,PETSC_VIEWER_STDOUT_(gcomm)); CHKERRQ(info);
	info = PetscPrintf(comm,"Hand-coded gradient\n");CHKERRQ(info);
	info = PetscObjectGetComm((PetscObject)g1,&gcomm);CHKERRQ(info);
	info = VecView(g1,PETSC_VIEWER_STDOUT_(gcomm));CHKERRQ(info);
	info = PetscPrintf(comm,"\n");CHKERRQ(info);
      }
      
      info = VecAXPY(g2,-1.0,g1); CHKERRQ(info); // g2 := g2-g1
      info = VecNorm(g1,NORM_2,&hcnorm); CHKERRQ(info);
      info = VecNorm(g2,NORM_2,&diffnorm); CHKERRQ(info);

      info = PetscPrintf(comm,"ratio ||fd-hc||/||hc|| = %G, difference ||fd-hc|| = %G\n",diffnorm/hcnorm,diffnorm);CHKERRQ(info);
      info = PetscPrintf(comm,"\n");CHKERRQ(info);
    }
  }
  if (fd->check_hessian) {
    info = TaoAppGetHessianMat(app,&A, &Apre); CHKERRQ(info);
    if (A != Apre) {
      SETERRQ(1,"Cannot test with alternative preconditioner");
    }
    info = PetscPrintf(comm,"Testing hand-coded hessian, if the ratio ||fd - hc|| / ||hc|| is\n"); CHKERRQ(info);
    info = PetscPrintf(comm,"0 (1.e-8), the hand-coded hessian is probably correct.\n"); CHKERRQ(info);

    if (!fd->complete_print) {
      info = PetscPrintf(comm,"Run with -tao_test_display to show difference\n");CHKERRQ(info);
      info = PetscPrintf(comm,"between hand-coded and finite difference hessian.\n");CHKERRQ(info);
    }
    for (i=0;i<3;i++) {
      if (i==1) {
	info = VecSet(x,-1.0); CHKERRQ(info);
      } else if (i==2) {
	info = VecSet(x,1.0); CHKERRQ(info);
      }
      info = TaoAppComputeHessian(app,x,&A,&A,&flg); CHKERRQ(info);
      if (!i) {
	info = MatConvert(A,MATSAME,MAT_INITIAL_MATRIX,&B); CHKERRQ(info);
      }
      info = TaoAppDefaultComputeHessian(app, x, &B, &B, &flg, PETSC_NULL); CHKERRQ(info);
      if (fd->complete_print) {
	info = PetscObjectGetComm((PetscObject)A,&gcomm); CHKERRQ(info);
	info = PetscPrintf(comm,"Finite difference hessian\n"); CHKERRQ(info);
	info = MatView(B,PETSC_VIEWER_STDOUT_(gcomm)); CHKERRQ(info);
	info = PetscPrintf(comm,"Hand-coded hessian\n"); CHKERRQ(info);
	info = MatView(A,PETSC_VIEWER_STDOUT_(gcomm)); CHKERRQ(info);
	info = PetscPrintf(comm,"\n"); CHKERRQ(info);
      }
      info = MatAXPY(B,-1.0,A,DIFFERENT_NONZERO_PATTERN); CHKERRQ(info);
      info = MatNorm(B,NORM_FROBENIUS,&diffnorm); CHKERRQ(info);
      info = MatNorm(A,NORM_FROBENIUS,&hcnorm); CHKERRQ(info);

      info = PetscPrintf(comm,"ratio ||fd-hc||/||hc|| = %G, difference ||fd-hc|| = %G\n",diffnorm/hcnorm,diffnorm);CHKERRQ(info);
      info = PetscPrintf(comm,"\n");CHKERRQ(info);
    }
    info = MatDestroy(B); CHKERRQ(info);
  }	
  

  TaoFunctionReturn(0);
  //  PetscFunctionReturn(PETSC_ERR_ARG_WRONGSTATE);

}

#undef __FUNCT__
#define __FUNCT__ "TaoSetUp_FD"
static int TaoSetUp_FD(TAO_SOLVER tao, void *solver)
{
  //TAO_FD *fdctx = (TAO_FD*)solver;
  TaoVec *X, *G;
  int info;
  TaoFunctionBegin;
  info = TaoGetSolution(tao,&X); CHKERRQ(info);
  info = X->Clone(&G); CHKERRQ(info);
  info = TaoSetLagrangianGradientVector(tao,G); CHKERRQ(info);
  info = TaoCheckFG(tao); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoDestroy_FD"
static int TaoDestroy_FD(TAO_SOLVER tao, void *solver)
{
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoSetOptions_FD"
static int TaoSetOptions_FD(TAO_SOLVER tao, void *solver) 
{
  TAO_FD *fd = (TAO_FD*)solver;
  int info;
  TaoFunctionBegin;
  info = PetscOptionsHead("Hand-coded gradient/hessian tester options");CHKERRQ(info);
  info = PetscOptionsName("-tao_test_display","Display difference between approximate and handcoded hessian","None",&fd->complete_print);CHKERRQ(info);
  info = PetscOptionsName("-tao_test_gradient","Test Hand-coded gradient against finite difference gradient","None",&fd->check_gradient);CHKERRQ(info);
  info = PetscOptionsName("-tao_test_hessian","Test Hand-coded hessian against finite difference hessian","None",&fd->check_hessian);CHKERRQ(info);
  if (fd->check_gradient == PETSC_FALSE && fd->check_hessian == PETSC_FALSE) {
    fd->check_gradient = PETSC_TRUE;
  }
    
  info = PetscOptionsTail();CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "TaoCreate_FD"
int TaoCreate_FD(TAO_SOLVER tao)
{
  TAO_FD *fd;
  int info;
  TaoFunctionBegin;
  info = TaoNew(TAO_FD,&fd); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_FD)); CHKERRQ(info);
  info = TaoSetTaoSolveRoutine(tao,TaoSolve_FD, (void*)fd); CHKERRQ(info);
  info = TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_FD, TaoDestroy_FD); CHKERRQ(info);
  info = TaoSetTaoOptionsRoutine(tao, TaoSetOptions_FD); CHKERRQ(info);
  //  info = TaoSetTaoViewRoutine(tao,TaoView_FD); CHKERRQ(info);
  fd->check_gradient=PETSC_TRUE;
  fd->check_hessian=PETSC_FALSE;
  fd->complete_print=PETSC_FALSE;
  TaoFunctionReturn(0);
}
EXTERN_C_END


#else
#undef __FUNCT__
#define __FUNCT__ "TaoCreate_FD"
EXTERN_C_BEGIN
int TaoCreate_FD(TAO_SOLVER tao)
{
  TaoFunctionBegin;
  SETERRQ(1,"tao_fd_test solver requires TAO_APPLICATION");
}
EXTERN_C_END

#endif
