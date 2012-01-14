#include "multigridmat.h"

#include "src/petsctao/vector/taovec_petsc.h"
#include "src/petsctao/indexset/taois_petsc.h"

TaoMultiGridMatPetsc::TaoMultiGridMatPetsc(Mat MM):TaoMatPetsc(MM){
  int i,nn=PETSCDAAPPMAXGRIDS;
  PetscFunctionBegin;
  for (i=0;i<nn;i++){
    this->grid[i].da=0;
    this->grid[i].H=0;
    this->grid[i].R=0;
    this->grid[i].RHS=0;
    this->grid[i].Interpolate=0;
    this->grid[i].CScale=0;
    this->grid[i].W3=0;
    this->grid[i].coloring=0;
    this->grid[i].mgrid=0;
  }
  this->nda=0;
  return;
}
 
TaoMultiGridMatPetsc::~TaoMultiGridMatPetsc(){
  this->TakeDown();
  return;
}
 

#undef __FUNCT__  
#define __FUNCT__ "TaoMultiGridMatPetsc::TakeDown"
int TaoMultiGridMatPetsc::TakeDown(){
  int i,info,nn=this->nda;
  PetscFunctionBegin;
  for (i=0;i<nn;i++){
    this->grid[i].da=0;
    this->grid[i].H=0;
    this->grid[i].R=0;
    this->grid[i].RHS=0;
    this->grid[i].Interpolate=0;
    this->grid[i].mgrid=0;
    this->grid[i].CScale=0;
    if (this->grid[i].W3){
      info=VecDestroy(this->grid[i].W3);CHKERRQ(info);
    }
    this->grid[i].W3=0;
  }
  this->nda=0;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoMultiGridMatPetsc::SetUpMultiGrid"
int TaoMultiGridMatPetsc::SetUpMultiGrid(GridCtx*dagrid,int nn){
  int i,info;
  PetscFunctionBegin;
  info=this->TakeDown(); CHKERRQ(info);
  this->nda=nn;
  this->ndamax=PETSCDAAPPMAXGRIDS;
  for (i=0;i<nn;i++){
    this->grid[i].H = dagrid[i].H;
    this->grid[i].R = dagrid[i].R;
    this->grid[i].RHS = dagrid[i].RHS;
    this->grid[i].Interpolate = dagrid[i].Interpolate;
    this->grid[i].CScale = dagrid[i].CScale;
    this->grid[i].W3=0;
    info=VecDuplicate(this->grid[i].R,&this->grid[i].W3); CHKERRQ(info);
  }
  info = this->SetMatrix(this->grid[nn-1].H,this->grid[nn-1].H,SAME_NONZERO_PATTERN); CHKERRQ(info);
  info=TaoSelectSubset(TaoMaskFullSpace); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoMultiGridMatPetsc::SetDiagonal"
int TaoMultiGridMatPetsc::SetDiagonal(TaoVec *tv)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);

  int i,info;
  Vec VCourse, VFine = pv->GetVec();

  PetscFunctionBegin;
  info = MatDiagonalSet(this->grid[this->nda-1].H,VFine,INSERT_VALUES); CHKERRQ(info);
  for (i=this->nda-1;i>0;i--){
    VCourse=this->grid[i-1].R;
    info=MatMultTranspose(this->grid[i].Interpolate,VFine,VCourse); CHKERRQ(info);
    info=VecPointwiseMult(VCourse,this->grid[i].CScale,VCourse); CHKERRQ(info);
    info=MatDiagonalSet(this->grid[i-1].H,VCourse,INSERT_VALUES); CHKERRQ(info);
    VFine=VCourse;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMultiGridMatPetsc::AddDiagonal"
int TaoMultiGridMatPetsc::AddDiagonal(TaoVec*tv) 
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);

  int i,info;
  Vec VCourse, VFine = pv->GetVec();

  PetscFunctionBegin;
  info = MatDiagonalSet(this->pm,VFine,ADD_VALUES); CHKERRQ(info);
  for (i=this->nda-1;i>0;i--){
    VCourse=this->grid[i-1].R;
    info=MatMultTranspose(this->grid[i].Interpolate,VFine,VCourse); CHKERRQ(info);
    info=VecPointwiseMult(VCourse,this->grid[i].CScale,VCourse); CHKERRQ(info);
    info=MatDiagonalSet(this->grid[i-1].H,VCourse,ADD_VALUES); CHKERRQ(info);
    VFine=VCourse;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoMultiGridMatPetsc::ShiftDiagonal"
int TaoMultiGridMatPetsc::ShiftDiagonal(double c){
  int i,info;
  PetscScalar cc=c;
  PetscFunctionBegin;
  info = MatShift(pm,cc); CHKERRQ(info);
  for (i=this->nda-1;i>0;i--){
    info = MatShift(this->grid[i-1].H,cc); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMultiGridMatPetsc::RowScale"
int TaoMultiGridMatPetsc::RowScale(TaoVec*tv)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);

  int i,info;
  Vec VCourse, VFine = pv->GetVec();

  PetscFunctionBegin;
  info = MatDiagonalScale(this->pm,VFine,PETSC_NULL); CHKERRQ(info);
  for (i=this->nda-1;i>0;i--){
    VCourse=this->grid[i-1].R;
    info=MatMultTranspose(this->grid[i].Interpolate,VFine,VCourse); CHKERRQ(info);
    info=VecPointwiseMult(VCourse,this->grid[i].CScale,VCourse); CHKERRQ(info);
    info=MatDiagonalScale(this->grid[i-1].H,VCourse,PETSC_NULL); CHKERRQ(info);
    VFine=VCourse;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMultiGridMatPetsc::ColScale"
int TaoMultiGridMatPetsc::ColScale(TaoVec*tv)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);

  int i,info;
  Vec VCourse,VFine = pv->GetVec();

  PetscFunctionBegin;
  info = MatDiagonalScale(this->pm,PETSC_NULL,VFine); CHKERRQ(info);
  for (i=this->nda-1;i>0;i--){
    VCourse=this->grid[i-1].R;
    info=MatMultTranspose(this->grid[i].Interpolate,VFine,VCourse); CHKERRQ(info);
    info=VecPointwiseMult(VCourse,this->grid[i].CScale,VCourse); CHKERRQ(info);
    info=MatDiagonalScale(this->grid[i-1].H,PETSC_NULL,VCourse); CHKERRQ(info);
    VFine=VCourse;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMultiGridMatPetsc::CreateReducedMatrix"
int TaoMultiGridMatPetsc::CreateReducedMatrix(TaoIndexSet*S1,TaoIndexSet*S2,TaoMat**MM){
  int info;
  Mat B,BPre;
  MatStructure flag;
  TaoMultiGridMatPetsc *MMM;

  PetscFunctionBegin;
  info = this->GetMatrix(&B,&BPre,&flag); CHKERRQ(info);
  MMM = new TaoMultiGridMatPetsc(B);
  info=MMM->SetUpMultiGrid(this->grid,this->nda); CHKERRQ(info);
  info = MMM->SetMatrix(B,BPre,flag); CHKERRQ(info);
  /*
    info = MatDuplicate(B,BB); CHKERRQ(info);
    info = MMM->SetMatrix(BB,BBPre,flag); CHKERRQ(info);
  */
  *MM=MMM;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMultiGridMatPetsc::SetReducedMatrix"
int TaoMultiGridMatPetsc::SetReducedMatrix(TaoMat*M,TaoIndexSet*S1,TaoIndexSet*S2){

  int i,info;
  TaoIndexSetPetsc *TRows=((TaoIndexSetPetsc *)S1);
  TaoIndexSetPetsc *TCols=((TaoIndexSetPetsc *)S2);
  TaoMatPetsc *MM=(TaoMatPetsc *)M;
  Vec CFine,CCourse,RFine,RCourse;
  Mat A,Apre,B,Bpre;

  PetscFunctionBegin;
  info=MM->GetMatrix(&A,&Apre,0); CHKERRQ(info);
  info=this->GetMatrix(&B,&Bpre,0); CHKERRQ(info);
  info=TRows->GetMask(&RFine); CHKERRQ(info);
  info=TCols->GetMask(&CFine); CHKERRQ(info);
  if (A!=B){
    info=MatCopy(A,B,SAME_NONZERO_PATTERN); CHKERRQ(info);
  }
  info=MatGetDiagonal(B,this->grid[this->nda-1].W3); CHKERRQ(info);
  info=MatDiagonalScale(B,RFine,CFine); CHKERRQ(info);
  info=MatDiagonalSet(B,this->grid[this->nda-1].W3,INSERT_VALUES); CHKERRQ(info);
  if (B!=Bpre){
    info=MatCopy(Apre,Bpre,SAME_NONZERO_PATTERN); CHKERRQ(info);
    info=MatGetDiagonal(Bpre,this->grid[this->nda-1].W3); CHKERRQ(info);
    info=MatDiagonalScale(Bpre,RFine,CFine); CHKERRQ(info);
    info=MatDiagonalSet(Bpre,this->grid[this->nda-1].W3,INSERT_VALUES); CHKERRQ(info);
  }
  for (i=this->nda-1;i>0;i--){
    CCourse=this->grid[i-1].R;
    RCourse=this->grid[i-1].RHS;
    info=MatMultTranspose(this->grid[i].Interpolate,CFine,CCourse); CHKERRQ(info);
    info=MatMultTranspose(this->grid[i].Interpolate,RFine,RCourse); CHKERRQ(info);
    info=VecPointwiseMult(CCourse,this->grid[i].CScale,CCourse); CHKERRQ(info);
    info=VecPointwiseMult(RCourse,this->grid[i].CScale,RCourse); CHKERRQ(info);
    info=MatGetDiagonal(this->grid[i-1].H,this->grid[i-1].W3); CHKERRQ(info);
    info=MatDiagonalScale(this->grid[i-1].H,CCourse,RCourse); CHKERRQ(info);
    info=MatDiagonalSet(this->grid[i-1].H,this->grid[i-1].W3,INSERT_VALUES); CHKERRQ(info);
    RFine=RCourse;
    CFine=CCourse;
  }
  PetscFunctionReturn(0);
}

