#include "tvecdouble.h"
#include "tao_general.h"
#include "stdio.h"

TaoVecDoubleArray::TaoVecDoubleArray(TaoInt nn):TaoVec(),n(nn){
  v=new double[nn];
  dallocated=1;
  return;
}

TaoVecDoubleArray::TaoVecDoubleArray(TaoInt nn, double *dd):TaoVec(),n(nn){
  v=dd;
  dallocated=0;
  return;
}

int TaoVecDoubleArray::Clone( TaoVec** tv ){

  *tv = new TaoVecDoubleArray(this->n);
  int info = (*tv)->CopyFrom(this);CHKERRQ(info);
  return 0;
}

int TaoVecDoubleArray::GetArray(TaoScalar **dptr, TaoInt *nn){
  if (sizeof(TaoScalar)==sizeof(double)){
    *dptr=(TaoScalar*)v;
    *nn=n;
  }
  else{
    return 1;
  }
  return 0;
}

int TaoVecDoubleArray::RestoreArray(TaoScalar **dptr, TaoInt *nn){
  *dptr=0;
  *nn=0;
  return 0;
}

int TaoVecDoubleArray::GetDoubles(double **dptr, TaoInt *nn){
  *dptr=v;
  *nn=n;
  return 0;
}

int TaoVecDoubleArray::RestoreDoubles(double **dptr, TaoInt *nn){
  *dptr=0;
  *nn=0;
  return 0;
}

int TaoVecDoubleArray::GetDimension(TaoInt *nn){
  *nn=n;
  return 0;
}

int TaoVecDoubleArray::Compatible(TaoVec *tv, TaoTruth *flag){
  TaoInt nn;
  int info;
  double *dptr;
  TaoVecDoubleArray* vv =  (TaoVecDoubleArray*)(tv);

  info = vv->GetData(&dptr,&nn);
  if (info==0 && nn == n) *flag=TAO_TRUE;
  else *flag=TAO_FALSE;
  return 0;
}

int TaoVecDoubleArray::View(){
  for (TaoInt i=0;i<n;++i)
    printf(" %4.2e \n ",v[i]);
  return 0;
}





