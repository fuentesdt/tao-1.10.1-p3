#include "tvecsingle.h"
#include "tao_general.h"
#include "stdio.h"

TaoVecFloatArray::TaoVecFloatArray(TaoInt nn):TaoVec(),n(nn){
  v=new float[nn];
  fallocated=1;
  return;
}

TaoVecFloatArray::TaoVecFloatArray(TaoInt nn, float *ff):TaoVec(),n(nn){
  v=ff;
  fallocated=0;
  return;
}

int TaoVecFloatArray::Clone( TaoVec** tv ){

  *tv = new TaoVecFloatArray(this->n);
  int info = (*tv)->CopyFrom(this);CHKERRQ(info);
  return 0;
}

int TaoVecFloatArray::GetFloats(float **fptr, TaoInt *nn){
  *fptr=v;
  *nn=n;
  return 0;
}

int TaoVecFloatArray::RestoreFloats(float **fptr, TaoInt *nn){
  *fptr=0;
  *nn=0;
  return 0;
}

int TaoVecFloatArray::GetDimension(TaoInt *nn){
  *nn=n;
  return 0;
}

int TaoVecFloatArray::Compatible(TaoVec *tv, TaoTruth *flag){
  TaoInt nn;
  int info;
  float *fptr;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);

  info = vv->GetData(&fptr,&nn);
  if (info==0 && nn == n) *flag=TAO_TRUE;
  else *flag=TAO_FALSE;
  return 0;
}

int TaoVecFloatArray::View(){
  for (TaoInt i=0;i<n;++i)
    printf(" %4.2e \n ",v[i]);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::SetToConstant"
int TaoVecFloatArray::SetToConstant( double c ){
  TaoInt i,nn;
  int info;
  float *tptr;

  TaoFunctionBegin;
  info = this->GetFloats(&tptr,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ tptr[i]=c; }
  info = this->RestoreFloats(&tptr,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::SetToZero"
int TaoVecFloatArray::SetToZero(){
  TaoInt i,nn;
  int info;
  float *tptr;

  TaoFunctionBegin;
  info = this->GetFloats(&tptr,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ tptr[i]=0; }
  info = this->RestoreFloats(&tptr,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::CopyFrom"
int TaoVecFloatArray::CopyFrom( TaoVec* tv ){
  TaoInt i,nn1,nn2;
  int info;
  float *tptr1,*tptr2;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&tptr2,&nn2);CHKERRQ(info);

  for (i=0;i<nn1;i++){ tptr1[i]=tptr2[i]; }

  info = vv->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::ScaleCopyFrom"
int TaoVecFloatArray::ScaleCopyFrom( double a, TaoVec* tv ){
  TaoInt i,nn1,nn2;
  int info;
  float *tptr1,*tptr2;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  for (i=0;i<nn1;i++){ tptr1[i]=a*tptr2[i]; }
  info = vv->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::NormInfinity"
int TaoVecFloatArray::NormInfinity(double *vnorm){
  TaoInt i,nn;
  int info;
  float dd=0, *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){
    if (vv[i]<0 && vv[i]<-dd) dd=-vv[i];
    else if (vv[i]>0 && vv[i]>dd) dd=vv[i];
  }
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  *vnorm=dd;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Norm1"
int TaoVecFloatArray::Norm1(double *vnorm){
  TaoInt i,nn;
  int info;
  float dd=0, *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (vv[i]<0) dd-=vv[i]; 
    else dd+=vv[i];
  }
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Norm2"
int TaoVecFloatArray::Norm2(double *vnorm){
  TaoInt i,nn;
  int info;
  float dd=0, *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) dd+=vv[i]*vv[i];
  *vnorm=sqrt(dd);
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Norm2squared"
int TaoVecFloatArray::Norm2squared(double *vnorm2){
  TaoInt i,nn;
  int info;
  float dd=0, *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) dd+=vv[i]*vv[i];
  *vnorm2=dd;
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Scale"
int TaoVecFloatArray::Scale( double alpha ){
  TaoInt i,nn;
  int info;
  float *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) vv[i]*=alpha;
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Axpy"
int TaoVecFloatArray::Axpy( double alpha, TaoVec* tv ){
  TaoInt i,nn1,nn2;
  int info;
  float *tptr1,*tptr2;
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tv);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = xx->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  for (i=0;i<nn1;i++){ tptr1[i]+= alpha * tptr2[i]; }
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = xx->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Axpby"
int TaoVecFloatArray::Axpby( double alpha, TaoVec* tv, double beta ){
  TaoInt i,nn1,nn2;
  int info;
  float *tptr1,*tptr2;
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tv);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = xx->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  for (i=0;i<nn1;i++){ tptr1[i] = beta * tptr1[i] + alpha * tptr2[i]; }
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = xx->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Aypx"
int TaoVecFloatArray::Aypx( double alpha, TaoVec* tv ){
  TaoInt i,nn1,nn2;
  int info;
  float *tptr1,*tptr2;
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tv);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = xx->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  for (i=0;i<nn1;i++){ tptr1[i] = alpha * tptr1[i] + tptr2[i]; }
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = xx->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  TaoFunctionReturn(0);
}
 
#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::AddConstant"
int TaoVecFloatArray::AddConstant( double alpha ){
  TaoInt i,nn;
  int info;
  float *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) vv[i]+=alpha;
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Dot"
int TaoVecFloatArray::Dot( TaoVec* tv, double *vDotv ){
  TaoInt i,nn1,nn2;
  int info;
  float dd=0,*tptr1,*tptr2;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  for (i=0;i<nn1;i++) dd+=tptr1[i]*tptr2[i];
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  *vDotv=dd;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Negate"
int TaoVecFloatArray::Negate(){ 
  TaoInt i,nn;
  int info;
  float *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) vv[i]=-vv[i];
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Reciprocal"
int TaoVecFloatArray::Reciprocal(){ 
  TaoInt i,nn;
  int info;
  float *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (vv[i]!=0)  vv[i]= 1.0/vv[i];
    else vv[i]=TAO_INFINITY;
  }
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Sqrt"
int TaoVecFloatArray::Sqrt(){ 
  TaoInt i,nn;
  int info;
  float *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (vv[i] >= 0)  vv[i]= sqrt(vv[i]);
    else vv[i]=TAO_INFINITY;
  }
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Pow"
int TaoVecFloatArray::Pow(double p){ 
  TaoInt i,nn;
  int info;
  float *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (vv[i] >= 0)  vv[i]= pow((float)vv[i], (float)p);
    else vv[i]=TAO_INFINITY;
  }
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::PointwiseMultiply"
int TaoVecFloatArray::PointwiseMultiply( TaoVec* tv, TaoVec* tw ){
  TaoInt i,nn1,nn2,nn3;
  int info;
  float *tptr1,*tptr2,*tptr3;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);
  TaoVecFloatArray* ww =  (TaoVecFloatArray*)(tw);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->GetFloats(&tptr3,&nn3);CHKERRQ(info);
  for (i=0;i<nn1;i++) tptr1[i]=tptr2[i] * tptr3[i];
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->RestoreFloats(&tptr3,&nn3);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::PointwiseDivide"
int TaoVecFloatArray::PointwiseDivide( TaoVec* tv , TaoVec* tw){
  TaoInt i,nn1,nn2,nn3;
  int info;
  float *tptr1,*tptr2,*tptr3;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);
  TaoVecFloatArray* ww =  (TaoVecFloatArray*)(tw);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->GetFloats(&tptr3,&nn3);CHKERRQ(info);

  for (i=0;i<nn1;i++) tptr1[i]=tptr2[i] / tptr3[i];
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->RestoreFloats(&tptr3,&nn3);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Median"
int TaoVecFloatArray::Median( TaoVec* tv, TaoVec* tw, TaoVec* tx){
  TaoInt i,nn1,nn2,nn3,nn4;
  int info;
  float *tptr1,*tptr2,*tptr3,*tptr4;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);
  TaoVecFloatArray* ww =  (TaoVecFloatArray*)(tw);
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tx);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->GetFloats(&tptr3,&nn3);CHKERRQ(info);
  info = xx->GetFloats(&tptr4,&nn4);CHKERRQ(info);

  for (i=0;i<nn1;i++){
    if (tptr2[i]<=tptr3[i] && tptr3[i] <= tptr4[i]){
      tptr1[i]=tptr3[i];
    } else if (tptr4[i]<=tptr3[i] && tptr3[i] <= tptr2[i]){
      tptr1[i]=tptr3[i];
    } else if (tptr3[i]<=tptr2[i] && tptr2[i] <= tptr4[i]){
      tptr1[i]=tptr2[i];
    } else if (tptr4[i]<=tptr2[i] && tptr2[i] <= tptr3[i]){
      tptr1[i]=tptr2[i];
    } else {
      tptr1[i]=tptr4[i];
    }
  }
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->RestoreFloats(&tptr3,&nn3);CHKERRQ(info);
  info = xx->RestoreFloats(&tptr4,&nn4);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::PointwiseMinimum"
int TaoVecFloatArray::PointwiseMinimum( TaoVec* tv, TaoVec* tw){
  TaoInt i,nn1,nn2,nn3;
  int info;
  float *tptr1,*tptr2,*tptr3;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);
  TaoVecFloatArray* ww =  (TaoVecFloatArray*)(tw);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->GetFloats(&tptr3,&nn3);CHKERRQ(info);

  for (i=0;i<nn1;i++) tptr1[i] = TaoMin( tptr2[i] , tptr3[i]);

  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->RestoreFloats(&tptr3,&nn3);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::PointwiseMaximum"
int TaoVecFloatArray::PointwiseMaximum( TaoVec* tv, TaoVec* tw){
  TaoInt i,nn1,nn2,nn3;
  int info;
  float *tptr1,*tptr2,*tptr3;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);
  TaoVecFloatArray* ww =  (TaoVecFloatArray*)(tw);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->GetFloats(&tptr3,&nn3);CHKERRQ(info);

  for (i=0;i<nn1;i++) tptr1[i] = TaoMax( tptr2[i] , tptr3[i]);

  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = vv->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  info = ww->RestoreFloats(&tptr3,&nn3);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Waxpby"
int TaoVecFloatArray::Waxpby  ( double a, TaoVec* tx, double b, TaoVec* ty){
  TaoInt i,nn1,nn2,nn3;
  int info;
  float *tptr1,*tptr2,*tptr3;
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tx);
  TaoVecFloatArray* yy =  (TaoVecFloatArray*)(ty);

  TaoFunctionBegin;
  info = this->GetFloats(&tptr1,&nn1);CHKERRQ(info);
  info = xx->GetFloats(&tptr2,&nn2);CHKERRQ(info);
  info = yy->GetFloats(&tptr3,&nn3);CHKERRQ(info);
  if (nn1!=nn2 || nn2!=nn3) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++){ tptr1[i] = a * tptr2[i] + b * tptr3[i]; }
  info = this->RestoreFloats(&tptr1,&nn1);CHKERRQ(info);
  info = xx->RestoreFloats(&tptr2,&nn2);CHKERRQ(info);
  info = yy->RestoreFloats(&tptr3,&nn3);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::AbsoluteValue"
int TaoVecFloatArray::AbsoluteValue(){
  TaoInt i,nn;
  int info;
  float *vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    vv[i]= TaoAbsScalar(vv[i]);
  }
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::MinElement"
int TaoVecFloatArray::MinElement(double*val){
  TaoInt i,nn;
  int info;
  float dd=TAO_INFINITY,*vv;

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (vv[i]<dd) dd=vv[i];
  }
  info = this->RestoreFloats(&vv,&nn);CHKERRQ(info);
  *val = dd;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::StepMax"
int TaoVecFloatArray::StepMax(TaoVec*tv,double*step1){
  TaoInt i,nn1,nn2;
  int info;
  float *xx,*dx;
  float stepmax1=TAO_INFINITY;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);

  TaoFunctionBegin;
  info = this->GetFloats(&xx,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&dx,&nn2);CHKERRQ(info);
  
  for (i=0;i<nn1;i++){
    if (xx[i] < 0){
      TaoFunctionReturn(1);
    } else if (dx[i]<0){ 
      stepmax1=TaoMin(stepmax1,-xx[i]/dx[i]);
    }
  }
  *step1=stepmax1;

  info = vv->RestoreFloats(&dx,&nn2);CHKERRQ(info);
  info = this->RestoreFloats(&xx,&nn1);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::StepMax2"
int TaoVecFloatArray::StepMax2(TaoVec*tv,TaoVec*txl,TaoVec*txu, double*step2){
  TaoInt i,nn1,nn2;
  int info;
  float *xx,*dx,*xl,*xu;
  float stepmax2=0;
  TaoVecFloatArray* vv =  (TaoVecFloatArray*)(tv);
  TaoVecFloatArray* xxll =  (TaoVecFloatArray*)(txl);
  TaoVecFloatArray* xxuu =  (TaoVecFloatArray*)(txu);

  TaoFunctionBegin;
  info = this->GetFloats(&xx,&nn1);CHKERRQ(info);
  info = vv->GetFloats(&dx,&nn2);CHKERRQ(info);
  info = xxll->GetFloats(&xl,&nn2);CHKERRQ(info);
  info = xxuu->GetFloats(&xu,&nn2);CHKERRQ(info);
  
  for (i=0;i<nn1;i++){
    if (dx[i] > 0){
      stepmax2=TaoMax(stepmax2,(xu[i]-xx[i])/dx[i]);      
    } else if (dx[i]<0){ 
      stepmax2=TaoMax(stepmax2,(xl[i]-xx[i])/dx[i]);
    }
  }
  info = this->RestoreFloats(&xx,&nn1);CHKERRQ(info);
  info = vv->RestoreFloats(&dx,&nn2);CHKERRQ(info);
  info = xxll->RestoreFloats(&xl,&nn2);CHKERRQ(info);
  info = xxuu->RestoreFloats(&xu,&nn2);CHKERRQ(info);
  
  *step2=stepmax2;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::BoundGradientProjection"
int TaoVecFloatArray::BoundGradientProjection(TaoVec*tg,TaoVec*txl,TaoVec*tx, TaoVec*txu){
  TaoInt i,nn1,nn2,nn3,nn4,nn5;
  int info;
  float *xptr,*xlptr,*xuptr,*gptr,*gpptr;
  TaoVecFloatArray* gg =  (TaoVecFloatArray*)(tg);
  TaoVecFloatArray* xxll =  (TaoVecFloatArray*)(txl);
  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tx);
  TaoVecFloatArray* xxuu =  (TaoVecFloatArray*)(txu);

  TaoFunctionBegin;
  info = this->GetFloats(&gpptr,&nn1);CHKERRQ(info);
  info = gg->GetFloats(&gptr,&nn2);CHKERRQ(info);
  info = xxll->GetFloats(&xlptr,&nn3);CHKERRQ(info);
  info = xx->GetFloats(&xptr,&nn4);CHKERRQ(info);
  info = xxuu->GetFloats(&xuptr,&nn5);CHKERRQ(info);
  if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4 || nn4!=nn5) {TaoFunctionReturn(1);}

  for (i=0; i<nn1; i++){

    gpptr[i] = gptr[i];
    if (gpptr[i]>0 && xptr[i]<=xlptr[i]){
      gpptr[i] = 0;
    } else if (gpptr[i]<0 && xptr[i]>=xuptr[i]){
      gpptr[i] = 0;
    }
  }
  info = this->RestoreFloats(&gpptr,&nn1);CHKERRQ(info);
  info = gg->RestoreFloats(&gptr,&nn2);CHKERRQ(info);
  info = xxll->RestoreFloats(&xlptr,&nn3);CHKERRQ(info);
  info = xx->RestoreFloats(&xptr,&nn4);CHKERRQ(info);
  info = xxuu->RestoreFloats(&xuptr,&nn5);CHKERRQ(info);

  TaoFunctionReturn(0);
}

inline static float fischer(float a, float b) 
{

#ifdef TODD

  if (TaoAbsScalar(a) > TaoAbsScalar(b)) {
    return sqrt(a*a + b*b) - a - b;
  }
  return sqrt(a*a + b*b) - b - a;

#else

   // Method suggested by Bob Vanderbei

   if (a + b <= 0) {
     return sqrt(a*a + b*b) - (a + b);
   }
   return -2.0*a*b / (sqrt(a*a + b*b) + (a + b));

#endif

}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::Fischer"
int TaoVecFloatArray::Fischer(TaoVec *tx, TaoVec *tf, TaoVec *tl, TaoVec *tu)
{
  TaoInt i,nn1,nn2,nn3,nn4,nn5;
  int info;
  float *vv,*x,*f,*l,*u;

  TaoVecFloatArray* xx =  (TaoVecFloatArray*)(tx);
  TaoVecFloatArray* ff =  (TaoVecFloatArray*)(tf);
  TaoVecFloatArray* ll =  (TaoVecFloatArray*)(tl);
  TaoVecFloatArray* uu =  (TaoVecFloatArray*)(tu);

  TaoFunctionBegin;
  info = this->GetFloats(&vv,&nn1);CHKERRQ(info);
  info = xx->GetFloats(&x,&nn2);CHKERRQ(info);
  info = ff->GetFloats(&f,&nn3);CHKERRQ(info);
  info = ll->GetFloats(&l,&nn4);CHKERRQ(info);
  info = uu->GetFloats(&u,&nn5);CHKERRQ(info);
  if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4 || nn4!=nn5) {
    TaoFunctionReturn(1);
  }

  for (i=0;i<nn1;i++) {

    if ((l[i] <= -TAO_INFINITY) && (u[i] >= TAO_INFINITY)) {
      vv[i] = -f[i];
    } 
    else if (l[i] <= -TAO_INFINITY) {
      vv[i] = -fischer(u[i] - x[i], -f[i]);
    } 
    else if (u[i] >=  TAO_INFINITY) {
      vv[i] =  fischer(x[i] - l[i],  f[i]);
    } 
    else if (l[i] == u[i]) {
      vv[i] = l[i] - x[i];
    }
    else {
      vv[i] =  fischer(u[i] - x[i], -f[i]);
      vv[i] =  fischer(x[i] - l[i],  vv[i]);
    }

  }

  info = this->RestoreFloats(&vv,&nn1);CHKERRQ(info);
  info = xx->RestoreFloats(&x,&nn2);CHKERRQ(info);
  info = ff->RestoreFloats(&f,&nn3);CHKERRQ(info);
  info = ll->RestoreFloats(&l,&nn4);CHKERRQ(info);
  info = uu->RestoreFloats(&u,&nn5);CHKERRQ(info);

  TaoFunctionReturn(0);
}

inline static float sfischer(float a, float b, float c)
{

#ifdef TODD

  if (TaoAbsScalar(a) > TaoAbsScalar(b)) {
    return sqrt(a*a + b*b + 2.0*c*c) - a - b;
  }
  return sqrt(a*a + b*b + 2.0*c*c) - b - a;

#else

   // Method suggested by Bob Vanderbei

   if (a + b <= 0) {
     return sqrt(a*a + b*b + 2.0*c*c) - (a + b);
   }
   return 2.0*(c*c - a*b) / (sqrt(a*a + b*b + 2.0*c*c) + (a + b));

#endif

}

#undef __FUNCT__
#define __FUNCT__ "TaoVecFloatArray::SFischer"
int TaoVecFloatArray::SFischer(TaoVec *tx, TaoVec *tf, TaoVec *tl, TaoVec *tu, double mu)
{
  TaoInt i, nn1, nn2, nn3, nn4, nn5;
  int info;
  float *vv, *x, *f, *l, *u;

  TaoVecFloatArray *xx = (TaoVecFloatArray *)(tx);
  TaoVecFloatArray *ff = (TaoVecFloatArray *)(tf);
  TaoVecFloatArray *ll = (TaoVecFloatArray *)(tl);
  TaoVecFloatArray *uu = (TaoVecFloatArray *)(tu);

  TaoFunctionBegin;

  if ((mu >= -TAO_EPSILON) && (mu <= TAO_EPSILON)) {
    Fischer(tx, tf, tl, tu);
  }
  else {
    info = this->GetFloats(&vv, &nn1); CHKERRQ(info);
    info = xx->GetFloats(&x, &nn2); CHKERRQ(info);
    info = ff->GetFloats(&f, &nn3); CHKERRQ(info);
    info = ll->GetFloats(&l, &nn4); CHKERRQ(info);
    info = uu->GetFloats(&u, &nn5); CHKERRQ(info);

    if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4 || nn4!=nn5) {
      TaoFunctionReturn(1);
    }

    for (i = 0; i < nn1; ++i) {
      if ((l[i] <= -TAO_INFINITY) && (u[i] >= TAO_INFINITY)) {
        vv[i] = -f[i] - mu*x[i];
      }
      else if (l[i] <= -TAO_INFINITY) {
        vv[i] = -sfischer(u[i] - x[i], -f[i], mu);
      }
      else if (u[i] >=  TAO_INFINITY) {
        vv[i] =  sfischer(x[i] - l[i],  f[i], mu);
      }
      else if (l[i] == u[i]) {
        vv[i] = l[i] - x[i];
      }
      else {
        vv[i] =  sfischer(u[i] - x[i], -f[i], mu);
        vv[i] =  sfischer(x[i] - l[i],  vv[i], mu);
      }
    }

    info = this->RestoreFloats(&vv, &nn1); CHKERRQ(info);
    info = xx->RestoreFloats(&x, &nn2); CHKERRQ(info);
    info = ff->RestoreFloats(&f, &nn3); CHKERRQ(info);
    info = ll->RestoreFloats(&l, &nn4); CHKERRQ(info);
    info = uu->RestoreFloats(&u, &nn5); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

