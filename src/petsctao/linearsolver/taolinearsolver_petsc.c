#include "tao_general.h"

#ifdef TAO_USE_PETSC

#include "taolinearsolver_petsc.h"

#include "../matrix/taomat_petsc.h"
#include "../vector/taovec_petsc.h"

TaoLinearSolverPetsc::TaoLinearSolverPetsc(KSP S):TaoLinearSolver(){
  linear_its=0;
  ksp=0;
  this->pkspviewer=PETSC_VIEWER_STDOUT_WORLD;
  this->SetKSP(S);
  return;
}

TaoLinearSolverPetsc::~TaoLinearSolverPetsc(){
  if (ksp){
    KSPDestroy(ksp);
  }
  return;
}


#undef __FUNCT__
#define __FUNCT__ "TaoWrapKSP"
/*@C
   TaoWrapKSP - Create a new TaoLinearSolver object using PETSc KSP.

   Input Parameter:
.  S -  a KSP

   Output Parameter:
.  SS - new TaoMat

   Note:  
   A TaoLinearSolverPetsc is an object with the methods of an abstract
   TaoLinearSolver object.  A TaoLinearSolverPetsc contains an implementation 
   of the TaoLinearSolver methods.  Routines using these vectors should 
   declare a pointer to a TaoLinearSolver, assign this pointer to the address 
   of a TaoLinearSolver object, use the pointer to invoke methods on the 
   object, and use this pointer as an argument when calling other routines.  
   This usage is different from the usage of a PETSc KSP.  In PETSc, 
   applications will typically declare a KSP, and pass it as an argument into 
   routines.  That is, applications will typically declare a pointer to a 
   TaoLinearSolver and use the pointer, or declare a KSP and use it directly.

   Level: developer

@*/
int TaoWrapKSP( KSP S, TaoLinearSolverPetsc ** SS){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(S,KSP_COOKIE,1);
  *SS = new TaoLinearSolverPetsc(S);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverGetKSP"
/*@C
   TaoLinearSolverGetKSP - If the TaoLinearSolver is of the TaoLinearSolverPetsc type, it gets the underlying PETSc KSP 

   Input Parameter:
+  TS - the TaoLinearSolver
-  S -  the address of KSP

   Output Parameter:
.  S -  address of the underlying PETSc KSP object


   Note:  
   This routine does not create a KSP.  It sets a pointer
   to the location of an existing KSP.

   Level: advanced

.seealso TaoVecGetPetscVec(), TaoMatGetPetscMat(), TaoAppGetKSP()
@*/
int TaoLinearSolverGetKSP( TaoLinearSolver *TS, KSP * S){
  PetscFunctionBegin;
  if (TS){
    *S=((TaoLinearSolverPetsc *)TS)->GetKSP();
    PetscValidHeaderSpecific(*S,KSP_COOKIE,2);
  } else {
    *S=0;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::SetKSP"
int TaoLinearSolverPetsc::SetKSP(KSP ksp2){
  int info;
  PetscFunctionBegin;
  if (ksp2){
    PetscValidHeaderSpecific(ksp2,KSP_COOKIE,1);
    PetscObjectReference((PetscObject)ksp2); 
  }
  if (ksp){
    info=KSPDestroy(ksp); CHKERRQ(info);
  }
  ksp=ksp2;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::GetKSP"
int TaoLinearSolverPetsc::GetKSP(KSP *ksp2){
  PetscFunctionBegin;
  if (ksp2){
    *ksp2=this->ksp;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::PreSolve"
int TaoLinearSolverPetsc::PreSolve(TaoMat* m1)
{
  TaoMatPetsc *MM = dynamic_cast <TaoMatPetsc *> (m1);
  Mat mm, mm_pre;
  MatStructure preflag;
  int info;

  PetscFunctionBegin;
  info = MM->GetMatrix(&mm, &mm_pre, &preflag); CHKERRQ(info);
  info = KSPSetOperators(ksp, mm, mm_pre, preflag); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::SetTrustRadius"
int TaoLinearSolverPetsc::SetTrustRadius(double rad)
{
  const KSPType ktype;
  int info;
  PetscTruth flg;

  PetscFunctionBegin;

  info = KSPGetType(ksp, &ktype); CHKERRQ(info);

  info = PetscStrcmp((char *)ktype, KSPNASH, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPNASHSetRadius(ksp, rad); CHKERRQ(info);
  }

  info = PetscStrcmp((char *)ktype, KSPSTCG, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPSTCGSetRadius(ksp, rad); CHKERRQ(info);
  }

  info = PetscStrcmp((char *)ktype, KSPGLTR, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPGLTRSetRadius(ksp, rad); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::GetNormDirection"
int TaoLinearSolverPetsc::GetNormDirection(double *norm_d)
{
  const KSPType ktype;
  int info;
  PetscTruth flg;

  PetscFunctionBegin;

  info = KSPGetType(ksp, &ktype); CHKERRQ(info);

  info = PetscStrcmp((char *)ktype, KSPNASH, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPNASHGetNormD(ksp, norm_d); CHKERRQ(info);
  }

  info = PetscStrcmp((char *)ktype, KSPSTCG, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPSTCGGetNormD(ksp, norm_d); CHKERRQ(info);
  }

  info = PetscStrcmp((char *)ktype, KSPGLTR, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPGLTRGetNormD(ksp, norm_d); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::GetObjFcn"
int TaoLinearSolverPetsc::GetObjFcn(double *o_fcn)
{
  const KSPType ktype;
  int info;
  PetscTruth flg;

  PetscFunctionBegin;

  info = KSPGetType(ksp, &ktype); CHKERRQ(info);

  info = PetscStrcmp((char *)ktype, KSPNASH, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPNASHGetObjFcn(ksp, o_fcn); CHKERRQ(info);
  }

  info = PetscStrcmp((char *)ktype, KSPSTCG, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPSTCGGetObjFcn(ksp, o_fcn); CHKERRQ(info);
  }

  info = PetscStrcmp((char *)ktype, KSPGLTR, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPGLTRGetObjFcn(ksp, o_fcn); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::GetMinEig"
int TaoLinearSolverPetsc::GetMinEig(double *e_min)
{
  const KSPType ktype;
  int info;
  PetscTruth flg;

  PetscFunctionBegin;

  *e_min = 0.0;

  info = KSPGetType(ksp, &ktype); CHKERRQ(info);
  info = PetscStrcmp((char *)ktype, KSPGLTR, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPGLTRGetMinEig(ksp, e_min); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::GetLambda"
int TaoLinearSolverPetsc::GetLambda(double *lambda)
{
  const KSPType ktype;
  int info;
  PetscTruth flg;

  PetscFunctionBegin;

  *lambda = 0.0;

  info = KSPGetType(ksp, &ktype); CHKERRQ(info);
  info = PetscStrcmp((char *)ktype, KSPGLTR, &flg); CHKERRQ(info);
  if (flg == PETSC_TRUE) { 	
    info = KSPGLTRGetLambda(ksp, lambda); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::Solve"
int TaoLinearSolverPetsc::Solve(TaoVec* tv, TaoVec* tw, TaoTruth *flag)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *pw = dynamic_cast <TaoVecPetsc *> (tw);
  int info;
  PetscInt its;

  PetscFunctionBegin;

  info = KSPSolve(ksp,pv->GetVec(),pw->GetVec()); CHKERRQ(info);
  info = KSPGetIterationNumber(ksp, &its); CHKERRQ(info);

  this->linear_its=PetscMax(its,-its);
  if (its>0) *flag=TAO_TRUE; 
  else *flag=TAO_FALSE;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::SolveTrustRegion"
int TaoLinearSolverPetsc::SolveTrustRegion(TaoVec *b, TaoVec *x, 
					   double r, TaoTruth *flg) 
{
  int info;

  PetscFunctionBegin;
  info = this->SetTrustRadius(r); CHKERRQ(info);
  info = this->Solve(b, x, flg); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::GetNumberIterations"
int TaoLinearSolverPetsc::GetNumberIterations(int * iters){
  PetscFunctionBegin;
  *iters=linear_its;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::SetTolerances"
int TaoLinearSolverPetsc::SetTolerances(double rtol, double atol,
					double dtol, int maxits) {
  int info;
  PetscFunctionBegin;
  info = KSPSetTolerances(ksp, rtol, atol, dtol, maxits); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::View"
int TaoLinearSolverPetsc::View(){
  int info;
  PetscFunctionBegin;
  info=KSPView(this->ksp,this->pkspviewer);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoLinearSolverPetsc::SetOptions"
int TaoLinearSolverPetsc::SetOptions(){
  int info;
  PetscFunctionBegin;
  info=KSPSetFromOptions(ksp); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#endif
