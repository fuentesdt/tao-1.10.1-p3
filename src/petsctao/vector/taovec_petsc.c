#include "tao_general.h"

#ifdef TAO_USE_PETSC

#include "taovec_petsc.h"
#include "../indexset/taois_petsc.h"

int VecMedian(Vec, Vec, Vec, Vec);
int VecPointwiseMax(Vec, Vec, Vec);
int VecPointwiseMin(Vec, Vec, Vec);
int VecFischer(Vec, Vec, Vec, Vec, Vec);
int VecSFischer(Vec, Vec, Vec, Vec, double, Vec);
int VecBoundProjection(Vec, Vec, Vec, Vec, Vec);
int VecCreateSubVec(Vec, IS, Vec *);
int VecISAXPY(Vec, double, Vec, IS);
int VecStepMax(Vec, Vec, double*);
int VecStepMax2(Vec, Vec,Vec, Vec, double*);
int VecCompare(Vec, Vec, PetscTruth *);
int VecBoundStepInfo(Vec, Vec, Vec, Vec, PetscReal *, PetscReal *, PetscReal *);
int VecPow(Vec, double);

static int petscminmaxevent=0;

#undef __FUNCT__
#define __FUNCT__ "TaoWrapPetscVec"
/*@C
   TaoWrapPetscVec - Creates a new TaoVec object using an existing
   PETSc Vec.

   Input Parameter:
+  V -  a PETSc vector
-  TV -  the address of a pointer to a TaoVecPetsc

   Output Parameter:
.  TV - pointer to a new TaoVecPetsc

   Note:  A TaoVecPetsc is an object with the methods of an abstract
   TaoVec object.  A TaoVecPetsc contains an implementation of the TaoVec
   methods.  Routines using these vectors should declare a pointer to 
   a TaoVec, assign this pointer to the address of a TaoVec object, 
   use the pointer to invoke methods on the object, and use this pointer
   as an argument when calling other routines.  This usage is different
   from the usage of a PETSc Vec.  In PETSc, applications will typically
   declare a Vec, and pass it as an argument into routines.  That is,
   applications will typically declare a pointer to a TaoVec and use the
   pointer, or declare a Vec and use it directly.
  
   Note:
   The user is responsible for destroying the Vec V, in addition to
   to TaoVecPetsc vector TV.  The Vec can be destroyed immediately
   after this routine.

   Level: developer

.seealso TaoVecGetPetscVec(), TaoVecDestroy()

@*/
int TaoWrapPetscVec( Vec V, TaoVecPetsc* *TV){
  int info;
  PetscFunctionBegin;
  *TV = new  TaoVecPetsc(V);
  info=(*TV)->SetVec(V); CHKERRQ(info);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "TaoVecGetPetscVec"                              
/*@C                                                               
   TaoVecGetPetscVec - Set the address of a Vec equal to the location
   of the underlying Vec structure in this TaoVecPetsc object.       
                                                                     
   Input Parameter:                                                  
+  TV - the TaoVecPetsc                                              
-  V -  the address of Vec                                           
                                                                     
   Output Parameter:                                                 
.  V -  the address of a specific Vec                                
                                                                     
   Note:                                                             
   This routine does not create a Vec.  It sets a pointer            
   to the location of an existing Vec.                               
                                                                     
                                                                     
   Note:                                                             
                                                                     
   VecObject is a void* pointer that can be cast to the underlying Vec object.
                                                                              
   Level: advanced                                                            
                                                                              
@*/                                                                           
int TaoVecGetPetscVec( TaoVec* TV, Vec *V){                                   
  PetscFunctionBegin;                                                         
  if (TV && V){         
    TaoVecPetsc *tvp = dynamic_cast<TaoVecPetsc*>(TV);
    if (tvp) {
      *V=(Vec)(tvp->GetVec());                                                  
    } else {
      *V=0;
    }
  } else if (V){                                                              
    *V=0;                                                                     
  }                                                                           
  PetscFunctionReturn(0);                                                     
}       

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::TaoVecPetsc"
TaoVecPetsc::TaoVecPetsc(Vec VV)
{
  pv=0;
  pvecviewer=0;
  if (VV) { 
    SetVec(VV);
  }
  if (petscminmaxevent==0){
    PetscLogEventRegister("PointwiseMinMax",VEC_COOKIE,&petscminmaxevent);
  }
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::~TaoVecPetsc"
TaoVecPetsc::~TaoVecPetsc()
{
  SetVec(0);
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::SetVec"
/*@C
   SetVec - Set or replace the underlying Vec object in this TaoVec.

   Input Parameter:
.  V -  a PETSc vector

   Note:
   The user us repsonsible for destroying the Vec V, in addition to
   to TaoVecPetsc vector TV.  The Vec can be destroyed immediately
   after this routine.

   Level: developer

.seealso TaoWrapPetscVec(), TaoVecDestroy()

@*/
int TaoVecPetsc::SetVec(Vec VV){
  int info;
  PetscFunctionBegin;
  if (VV){
    PetscValidHeaderSpecific(VV,VEC_COOKIE,1);
    PetscObjectReference((PetscObject)VV); 
    //    info = PetscInfo(VV,"Set the PETSc Vec in a TaoVec .\n"); CHKERRQ(info);
  }
  if (pv) {
    info=VecDestroy(pv); CHKERRQ(info);
  }
  pv=VV;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Clone"
int TaoVecPetsc::Clone(TaoVec**ntv){
  int info;
  TaoVecPetsc *nptv;
  Vec npv;
  PetscFunctionBegin;
  info=VecDuplicate(pv,&npv); CHKERRQ(info);
  info = TaoWrapPetscVec(npv,&nptv); CHKERRQ(info);
  *ntv=nptv;
  info = VecDestroy(npv); CHKERRQ(info);
  nptv->pvecviewer=pvecviewer;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Compatible"
int TaoVecPetsc::Compatible(TaoVec *tv, TaoTruth *flag){
  TaoVecPetsc *vv = (TaoVecPetsc *)(tv);
  PetscTruth flg=PETSC_TRUE;
  int info=0;

  PetscFunctionBegin;
  if (vv==0){
    *flag=TAO_FALSE;
    PetscFunctionReturn(0);
  }
  info = VecCompare(pv, vv->pv, &flg); CHKERRQ(info);
  if (flg==PETSC_FALSE){
    *flag=TAO_FALSE;
  } else {
    *flag=TAO_TRUE;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::SetToConstant"
int TaoVecPetsc::SetToConstant( double cc ){
  int info;
  PetscScalar c=cc;
  PetscFunctionBegin;
  info=VecSet(pv, c);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::SetToZero"
int TaoVecPetsc::SetToZero(){
  PetscScalar zero=0.0;
  int info;
  PetscFunctionBegin;
  info=VecSet(pv, zero);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::CopyFrom"
int TaoVecPetsc::CopyFrom(TaoVec* tv)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  int info;

  PetscFunctionBegin;
  info = VecCopy(mv->pv, pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::ScaleCopyFrom"
int TaoVecPetsc::ScaleCopyFrom(double a, TaoVec *tv)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  PetscScalar alpha=a, zero=0.0;
  int info;
  PetscFunctionBegin;
  info=VecAXPBY(pv,alpha,zero,mv->pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::NormInfinity"
int TaoVecPetsc::NormInfinity(double *vnorm){
  int info;
  PetscReal vv;
  PetscFunctionBegin;
  info=VecNorm(pv,NORM_INFINITY,&vv);  CHKERRQ(info);
  *vnorm=(double)vv;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Norm1"
int TaoVecPetsc::Norm1(double *vnorm){
  int info;
  PetscReal vv;
  PetscFunctionBegin;
  info=VecNorm(pv,NORM_1,&vv);  CHKERRQ(info);
  *vnorm=(double)vv;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Norm2"
int TaoVecPetsc::Norm2(double *vnorm){
  int info;
  PetscReal vv;
  PetscFunctionBegin;
  info=VecNorm(pv,NORM_2,&vv);  CHKERRQ(info);
  *vnorm=(double)vv;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Norm2squared"
int TaoVecPetsc::Norm2squared(double *vnorm){
  int info;
  PetscReal vv;
  PetscFunctionBegin;
  info=VecNorm(pv,NORM_2,&vv);  CHKERRQ(info);
  *vnorm=(double)(vv*vv);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Scale"
int TaoVecPetsc::Scale( double alpha ){
  int info;
  PetscScalar aalpha=alpha;
  PetscFunctionBegin;
  info=VecScale(pv, aalpha);  CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Axpy"
int TaoVecPetsc::Axpy(double alpha, TaoVec *tv)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  PetscScalar aalpha=alpha;
  int info;
  PetscFunctionBegin;
  info=VecAXPY(pv, aalpha, mv->pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Axpby"
int TaoVecPetsc::Axpby(double alpha, TaoVec *tv, double beta)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  PetscScalar aalpha = alpha, bbeta = beta;
  int info;
  PetscFunctionBegin;
  info=VecAXPBY(pv,aalpha,bbeta,mv->pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Aypx"
int TaoVecPetsc::Aypx(double alpha, TaoVec *tv)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  PetscScalar aalpha = alpha;
  int info;
  PetscFunctionBegin;
  info=VecAYPX(pv,aalpha,mv->pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}
 
#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::AddConstant"
int TaoVecPetsc::AddConstant( double cc )
{
  int info;
  PetscScalar c=cc;
  PetscFunctionBegin;
  info=VecShift(pv, c);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Dot"
int TaoVecPetsc::Dot( TaoVec* tv, double *vDotv )
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  PetscScalar c;
  int info;
  PetscFunctionBegin;
  info=VecDot(pv, mv->pv, &c);  CHKERRQ(info);
  *vDotv=c;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Negate"
int TaoVecPetsc::Negate(){ 
  PetscScalar m1=-1.0; 
  int info;
  PetscFunctionBegin;
  info=VecScale(pv, m1);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Reciprocal"
int TaoVecPetsc::Reciprocal(){ 
  int info;
  PetscFunctionBegin;
  info=VecReciprocal(pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Sqrt"
int TaoVecPetsc::Sqrt(){ 
  int info;
  PetscFunctionBegin;
  info=VecSqrt(pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Pow"
int TaoVecPetsc::Pow(double p){ 
  int info;
  PetscFunctionBegin;
  info=VecPow(pv, p);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::GetDimension"
int TaoVecPetsc::GetDimension(TaoInt *n){
  int info;
  PetscFunctionBegin;
  info=VecGetSize(pv,n);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::PointwiseMultiply"
int TaoVecPetsc::PointwiseMultiply(TaoVec *tv, TaoVec *tw)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *mw = dynamic_cast <TaoVecPetsc *> (tw);
  int info;
  PetscFunctionBegin;
  info=VecPointwiseMult(pv, mv->pv, mw->pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::PointwiseDivide"
int TaoVecPetsc::PointwiseDivide( TaoVec* tv , TaoVec* tw)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *mw = dynamic_cast <TaoVecPetsc *> (tw);
  int info;
  PetscFunctionBegin;
  info=VecPointwiseDivide(pv, mv->pv, mw->pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Median"
int TaoVecPetsc::Median( TaoVec* tv, TaoVec* tw, TaoVec* tx)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *mw = dynamic_cast <TaoVecPetsc *> (tw);
  TaoVecPetsc *mx = dynamic_cast <TaoVecPetsc *> (tx);
  int info;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscminmaxevent,0,0,0,0); CHKERRQ(info);
  info=VecMedian(mv->pv,mw->pv,mx->pv,pv); CHKERRQ(info);
  info=PetscLogEventEnd(petscminmaxevent,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::PointwiseMinimum"
int TaoVecPetsc::PointwiseMinimum( TaoVec* tv, TaoVec* tw)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *mw = dynamic_cast <TaoVecPetsc *> (tw);
  int info;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscminmaxevent,0,0,0,0); CHKERRQ(info);
  info=VecPointwiseMin(pv, mv->pv, mw->pv); CHKERRQ(info);
  info=PetscLogEventEnd(petscminmaxevent,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Fischer"
int TaoVecPetsc::Fischer(TaoVec *tx, TaoVec *tf, TaoVec *tl, TaoVec *tu) 
{
  TaoVecPetsc *xx = dynamic_cast <TaoVecPetsc *> (tx);
  TaoVecPetsc *ff = dynamic_cast <TaoVecPetsc *> (tf);
  TaoVecPetsc *ll = dynamic_cast <TaoVecPetsc *> (tl);
  TaoVecPetsc *uu = dynamic_cast <TaoVecPetsc *> (tu);
  int info;

  PetscFunctionBegin;
  info=VecFischer(xx->pv, ff->pv, ll->pv, uu->pv, pv); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::SFischer"
int TaoVecPetsc::SFischer(TaoVec *tx, TaoVec *tf,
                          TaoVec *tl, TaoVec *tu, double mu)
{
  TaoVecPetsc *xx = dynamic_cast <TaoVecPetsc *> (tx);
  TaoVecPetsc *ff = dynamic_cast <TaoVecPetsc *> (tf);
  TaoVecPetsc *ll = dynamic_cast <TaoVecPetsc *> (tl);
  TaoVecPetsc *uu = dynamic_cast <TaoVecPetsc *> (tu);
  int info;

  PetscFunctionBegin;
  if ((mu >= -TAO_EPSILON) && (mu <= TAO_EPSILON)) {
    Fischer(tx, tf, tl, tu);
  }
  else {
    info=VecSFischer(xx->pv, ff->pv, ll->pv, uu->pv, mu, pv); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::PointwiseMaximum"
int TaoVecPetsc::PointwiseMaximum( TaoVec* tv, TaoVec* tw)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *mw = dynamic_cast <TaoVecPetsc *> (tw);
  int info;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscminmaxevent,0,0,0,0); CHKERRQ(info);
  info=VecPointwiseMax(pv, mv->pv, mw->pv);  CHKERRQ(info);
  info=PetscLogEventEnd(petscminmaxevent,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0); 
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::Waxpby"
int TaoVecPetsc::Waxpby  ( double aa, TaoVec* tv, double bb, TaoVec* tw)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *mw = dynamic_cast <TaoVecPetsc *> (tw);
  int info;
  PetscScalar a=aa, b=bb, zero=0.0;

  PetscFunctionBegin;
  if (a==1){
    info = VecWAXPY(pv,b,mw->pv,mv->pv); CHKERRQ(info);
  } 
  else if (b==1) {
    info = VecWAXPY(pv,a,mv->pv,mw->pv); CHKERRQ(info);
  } 
  else if (a==0) {
    info = VecSet(pv, zero); CHKERRQ(info);
    info = VecAXPY(pv, b, mw->pv); CHKERRQ(info);
  } 
  else if (b==0) {
    info = VecSet(pv, zero); CHKERRQ(info);
    info = VecAXPY(pv, a, mv->pv); CHKERRQ(info);
  } 
  else {    
    info = VecCopy(mw->pv,pv); CHKERRQ(info);
    info = VecAXPBY(pv,a,b,mv->pv); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::AbsoluteValue"
int TaoVecPetsc::AbsoluteValue(){
  int info;
  PetscFunctionBegin;
  info=VecAbs(pv);  CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::min"
int TaoVecPetsc::MinElement(double*val){
  int info;
  PetscInt p;
  PetscScalar a;
  PetscFunctionBegin;
  info = VecMin(pv,&p, &a);CHKERRQ(info);
  *val=a;
  PetscFunctionReturn(0);  
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::GetArray"
int TaoVecPetsc::GetArray(TaoScalar **dptr, TaoInt *n){
  int info;
  PetscScalar *pptr;
  PetscFunctionBegin;
  if (sizeof(TaoScalar)==sizeof(PetscScalar)){
    info = VecGetLocalSize(pv,n);CHKERRQ(info);
    info = VecGetArray(pv,&pptr);CHKERRQ(info);
    *dptr=(TaoScalar*)pptr;
  } else {
    SETERRQ(1,"Incompatible data types");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::RestoreArray"
int TaoVecPetsc::RestoreArray(TaoScalar **dptr,TaoInt *n){
  int info;
  PetscScalar *pptr;
  PetscFunctionBegin;
  if (sizeof(TaoScalar)==sizeof(PetscScalar)){
    pptr=(PetscScalar*)(*dptr);
    info = VecRestoreArray(pv,&pptr); CHKERRQ(info);
    *dptr=PETSC_NULL;
    *n=0;
  } else {
    SETERRQ(1,"Incompatible data types");
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::SetReducedVec"
int TaoVecPetsc::SetReducedVec(TaoVec *tv, TaoIndexSet *ti)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoIndexSetPetsc *TSS = dynamic_cast <TaoIndexSetPetsc *> (ti);

  int info;
  PetscInt nn,nnn;
  TaoPetscISType type;

  PetscFunctionBegin;
  info = TSS->GetSize(&nn); CHKERRQ(info); 
  info = TSS->GetReducedType(&type); CHKERRQ(info);
  info = mv->GetDimension(&nnn); CHKERRQ(info); 
  if (type==TaoMaskFullSpace){
  } 
  else if (type==TaoMatrixFree){
    info = VecCopy(mv->pv,pv); CHKERRQ(info);
    info = VecISSetToConstant(TSS->GetIS(),0.0,pv);CHKERRQ(info);
  } 

  info=ReducedCopyFromFull(mv,TSS);CHKERRQ(info);

  PetscFunctionReturn(0);  
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::ReducedXPY"
int TaoVecPetsc::ReducedXPY(TaoVec *tv, TaoIndexSet *ti)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoIndexSetPetsc *TSS = dynamic_cast <TaoIndexSetPetsc *> (ti);

  VecScatter scatterit;
  int info;

  PetscFunctionBegin;

  info = TSS->GetReducedVecScatter(mv->pv,pv, &scatterit);CHKERRQ(info); 
  
  info = VecScatterBegin(scatterit,mv->pv,pv,ADD_VALUES,SCATTER_REVERSE);
  CHKERRQ(info); 
  info = VecScatterEnd(scatterit,mv->pv,pv,ADD_VALUES,SCATTER_REVERSE);
  CHKERRQ(info);

  PetscFunctionReturn(0);  
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::ReducedCopyFromFull"
int TaoVecPetsc::ReducedCopyFromFull(TaoVec *tv, TaoIndexSet *ti)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoIndexSetPetsc *TSS = dynamic_cast <TaoIndexSetPetsc *> (ti);

  int info;
  PetscScalar zero=0.0;
  VecScatter scatterit;

  PetscFunctionBegin;

  info = VecSet(pv, zero); CHKERRQ(info); 
  info = TSS->GetReducedVecScatter(mv->pv,pv, &scatterit);CHKERRQ(info); 
  info = VecScatterBegin(scatterit,mv->pv,pv,INSERT_VALUES,SCATTER_FORWARD);
  CHKERRQ(info); 
  info = VecScatterEnd(scatterit,mv->pv,pv,INSERT_VALUES,SCATTER_FORWARD);
  CHKERRQ(info);

  PetscFunctionReturn(0);  
}



#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::StepMax"
int TaoVecPetsc::StepMax(TaoVec *tv, double *step)
{
  TaoVecPetsc *mv = dynamic_cast <TaoVecPetsc *> (tv);
  PetscReal pstep;
  int info;
  PetscFunctionBegin;
  info = VecStepMax(pv, mv->pv, &pstep);CHKERRQ(info);
  *step=pstep;
  PetscFunctionReturn(0);  
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::StepBoundInfo"
int TaoVecPetsc::StepBoundInfo(TaoVec *txl ,TaoVec *txu, TaoVec *ts, 
                               double *bmin1, double *bmin2, double *bmax)
{
  TaoVecPetsc *mxl = dynamic_cast <TaoVecPetsc *> (txl);
  TaoVecPetsc *mxu = dynamic_cast <TaoVecPetsc *> (txu);
  TaoVecPetsc *ms = dynamic_cast <TaoVecPetsc *> (ts);

  int info;
  PetscReal p1,p2,p3;
  PetscFunctionBegin;
  info=VecBoundStepInfo(pv,mxl->pv,mxu->pv,ms->pv,&p1,&p2,&p3); CHKERRQ(info);
  *bmin1=p1; *bmin2=p2; *bmax=p3;
  PetscFunctionReturn(0);  
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::SetPetscViewer"
/*@C
   SetPetscViewer

   Input Parameter:
.  viewer - a viewer object to be used with View() and VecView()

   Level: advanced

@*/
int TaoVecPetsc::SetPetscViewer(PetscViewer viewer)
{
  PetscFunctionBegin;  
  pvecviewer= viewer;
  PetscFunctionReturn(0);  
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::View"
int TaoVecPetsc::View()
{
  int info;
  PetscFunctionBegin;  
  info=VecView(pv,pvecviewer);CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::BoundGradientProjection"
int TaoVecPetsc::BoundGradientProjection(TaoVec *tg, TaoVec *txl, 
				         TaoVec *tx, TaoVec *txu)
{
  TaoVecPetsc *mg = dynamic_cast <TaoVecPetsc *> (tg);
  TaoVecPetsc *mx = dynamic_cast <TaoVecPetsc *> (tx);
  TaoVecPetsc *mxl = dynamic_cast <TaoVecPetsc *> (txl);
  TaoVecPetsc *mxu = dynamic_cast <TaoVecPetsc *> (txu);
  int info;
  PetscFunctionBegin;
  info = PetscLogEventBegin(petscminmaxevent,0,0,0,0);CHKERRQ(info);
  info = VecBoundProjection(mg->pv,mx->pv,mxl->pv,mxu->pv,pv);CHKERRQ(info);
  info = PetscLogEventEnd(petscminmaxevent,0,0,0,0);CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecPetsc::CreateIndexSet"
int TaoVecPetsc::CreateIndexSet(TaoIndexSet**SSS)
{
  TaoIndexSetPetsc *SS;

  PetscFunctionBegin;  
  SS = new TaoIndexSetPetsc(pv);
  *SSS=SS;
  PetscFunctionReturn(0);  
}


#endif





