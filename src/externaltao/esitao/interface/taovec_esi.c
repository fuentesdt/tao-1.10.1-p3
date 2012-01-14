
#include "esi/ESI.h"
#include "tao_general.h"
#include "taovec_esi.h"
#include "stdio.h"

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::TaoVecESI"
TaoVecESI::TaoVecESI(esi::Vector<double,int> *esiv):TaoVec(),esivec(esiv){
  esiv->getLocalSize(this->nlocal); 
  esiv->addReference(); 
  this->VecObject=(void*)esiv; 
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::~TaoVecESI"
TaoVecESI::~TaoVecESI(){
  this->esivec->deleteReference();
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Compatible"
int TaoVecESI::Compatible(TaoVec* tv, TaoTruth *flag){
  int n1,n2;
  esi::ErrorCode info;
  esi::IndexSpace<int> *space1, *space2;
  TaoVecESI *that = (TaoVecESI *)(tv);
  info = this->esivec->getIndexSpace(space1); CHKERRQ(info);
  info = that->esivec->getIndexSpace(space2); CHKERRQ(info);

  info = space1->getGlobalSize(n1); CHKERRQ(info);
  info = space2->getGlobalSize(n2); CHKERRQ(info);
  if (n1!=n2){
    *flag= TAO_FALSE;
    TaoFunctionReturn(0);
  }
  info = space1->getLocalSize(n1); CHKERRQ(info);
  info = space2->getLocalSize(n2); CHKERRQ(info);
  if (n1!=n2){
    *flag= TAO_FALSE;
    TaoFunctionReturn(0);
  }
  /*
  info = space1->getLocalPartitionOffset(n1); if (info)  return TAO_FALSE;
  info = space2->getLocalPartitionOffset(n2); if (info)  return TAO_FALSE;
  if (n1!=n2) return TAO_FALSE;
  */
  *flag=TAO_TRUE;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Clone"
int TaoVecESI::Clone( TaoVec **ntv ){
  esi::ErrorCode info;
  esi::Vector<double,int> *nesi;

  TaoFunctionBegin;
  info=this->esivec->clone(nesi);  CHKERRQ(info);
  *ntv = new TaoVecESI(nesi);
  info=nesi->deleteReference(); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::GetDimension"
int TaoVecESI::GetDimension( int *n ){
  esi::ErrorCode info;

  TaoFunctionBegin;
  info=this->esivec->getGlobalSize(*n); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::SetToConstant"
int TaoVecESI::SetToConstant( double c ){
  double cc = c;
  esi::ErrorCode info;

  TaoFunctionBegin;
  info=this->esivec->put(cc); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::SetToZero"
int TaoVecESI::SetToZero(){
  esi::ErrorCode info;

  TaoFunctionBegin;
  info=esivec->put(0.0); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::CopyFrom"
int TaoVecESI::CopyFrom( TaoVec* tv ){
  TaoVecESI *vv = (TaoVecESI *)(tv);
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->copy(*vv->esivec); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::ScaleCopyFrom"
int TaoVecESI::ScaleCopyFrom( double a, TaoVec* tv ){
  TaoVecESI *vv = (TaoVecESI *)(tv);
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->copy(*vv->esivec); CHKERRQ(info);
  info=esivec->scale(a); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::NormInfinity"
int TaoVecESI::NormInfinity(double *vnorm){
  esi::ErrorCode info;
  double vvnorm;
  TaoFunctionBegin;
  info=esivec->normInfinity(vvnorm); CHKERRQ(info);
  *vnorm=vvnorm;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Norm1"
int TaoVecESI::Norm1(double *vnorm){
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->norm1(*vnorm); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Norm2"
int TaoVecESI::Norm2(double *vnorm){
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->norm2(*vnorm); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Norm2squared"
int TaoVecESI::Norm2squared(double *vnorm){
  esi::ErrorCode info;
  double vvnorm;
  TaoFunctionBegin;
  info=esivec->norm2squared(vvnorm); CHKERRQ(info);
  *vnorm=vvnorm;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Scale"
int TaoVecESI::Scale( TaoScalar alpha ){
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->scale(alpha); CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Axpy"
int TaoVecESI::Axpy( TaoScalar alpha, TaoVec* tv ){
  TaoVecESI *vv = (TaoVecESI *)(tv);
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->axpy(*vv->esivec,alpha);  CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Aypx"
int TaoVecESI::Aypx( TaoScalar alpha, TaoVec* x ){
  TaoVecESI *vv = (TaoVecESI *)(x);
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->aypx(alpha,*vv->esivec);  CHKERRQ(info);
  TaoFunctionReturn(0);
}
 

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Dot"
int TaoVecESI::Dot( TaoVec* tv, TaoScalar *vdotv ){
  TaoVecESI *vv = (TaoVecESI *)(tv);
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->dot(*vv->esivec,*vdotv); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Negate"
int TaoVecESI::Negate(){ 
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->scale(-1.0);  CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::PointwiseMultiply"
int TaoVecESI::PointwiseMultiply( TaoVec* tv, TaoVec* tw ){
  TaoVecESI *vv = (TaoVecESI *)(tv);
  TaoVecESI *ww = (TaoVecESI *)(tw);
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esivec->copy(*vv->esivec); CHKERRQ(info);
  info=esivec->scaleDiagonal(*ww->esivec); CHKERRQ(info);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::Waxpby"
int TaoVecESI::Waxpby  ( TaoScalar a, TaoVec* tv, TaoScalar b, TaoVec* tw){
  esi::ErrorCode info;
  TaoVecESI *vv = (TaoVecESI *)(tv);
  TaoVecESI *ww = (TaoVecESI *)(tw);

  TaoFunctionBegin;
  info=esivec->axpby(a,*vv->esivec,b,*ww->esivec); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::GetArray"
int TaoVecESI::GetArray(TaoScalar ** dptr, int * nn){
  esi::ErrorCode info;

  TaoFunctionBegin;
  info=esivec->getCoefPtrReadWriteLock(*dptr); CHKERRQ(info);
  *nn=this->nlocal;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::RestoreArray"
int TaoVecESI::RestoreArray(TaoScalar ** dptr, int * nn){
  esi::ErrorCode info;

  TaoFunctionBegin;
  info=esivec->releaseCoefPtrLock(*dptr); CHKERRQ(info);
  *dptr=0;
  TaoFunctionReturn(0);
}


#include <iostream.h>
#undef __FUNCT__
#define __FUNCT__ "TaoVecESI::View"
int TaoVecESI::View(){
  esi::ErrorCode info;
  esi::IndexSpace<int>* esimappart = NULL;
  info = esivec->getIndexSpace(esimappart); CHKERRQ(info);

   //we need to know how many local equations there are, and what the first
   //local equation is. (Since we happen to know that the map was built with
   //an indexBase of 0, then firstLocalEqn == 'getLocalPartitionOffset'.)
   int localSize, firstLocalEqn=0, localRank=0;
   info = esimappart->getLocalSize(localSize); CHKERRQ(info);
   info = esimappart->getLocalPartitionOffset(firstLocalEqn); CHKERRQ(info);
   info = esimappart->getLocalPartitionRank(localRank); CHKERRQ(info);

   double* coefs;
   info = esivec->getCoefPtrReadLock(coefs); CHKERRQ(info);

   for(int i=0; i<localSize; i++) {
      cout << localRank << ": "<<firstLocalEqn+i<<", "<<coefs[i]<<endl;
   }

   info = esivec->releaseCoefPtrLock(coefs); CHKERRQ(info);

   return(0);
}
