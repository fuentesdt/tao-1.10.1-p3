include "tao_general.h"
#include "taoesis.h"


#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::TaoESIIndexSet"
/* @C
   TaoESIIndexSet - Create a new TaoIndexSet object that works well with
   ESI objects.

   Input Parameter:
-  imax - the maximum local length of the index set (Should be equal to or greater than the local length of the vector)
+  SS - an Index Set

   Level: beginner

@ */
TaoESIIndexSet::TaoESIIndexSet(int imax):TaoIndexSet(){

  int size;

  iptr = new int[imax];
  nlocal=imax;

  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::TaoESIIndexSet"
TaoESIIndexSet::TaoESIIndexSet(int imax):TaoIndexSet(){

  int size;

  iptr = new int[imax];
  nlocal=imax;

  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::duplicate"
int TaoESIIndexSet::duplicate(TaoIndexSet**S){
  int info;
  IS is;

  TaoFunctionBegin;
  *S = new TaoESIIndexSet(nlocal, is);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::isSame"
int TaoESIIndexSet::isSame(TaoIndexSet*SS, TaoTruth*flg){
  int info;
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::whichEqual"
int TaoESIIndexSet::whichEqual(TaoVec*v1,TaoVec*v2){
  int i,info;
  int n,n2;
  double *vv1,*vv2;

  TaoFunctionBegin;
  info = v1->getArray(&vv1,&n); CHKERRQ(info);
  if (v1 == v2){
    vv2=vv1;
  } else {
    info = v2->getArray(&vv2,&n2); CHKERRQ(info);
  }

  if ( n != nlocalmax || n != n2 ){
    SETERRQ(1,"Vectors must be identically loaded over processors");
  }

  nlocal=0;
  for (i=0; i<n; i++){
    if (v1[i] == v2[i]) {same[i]=1; nlocal++}
    else same[i]=0;
  }
  
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::whichLessThan"
int TaoESIIndexSet::whichLessThan(TaoVec*v1,TaoVec*v2){
  int info;
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::whichGreaterThan"
int TaoESIIndexSet::whichGreaterThan(TaoVec*v1,TaoVec*v2){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::whichBetween"
int TaoESIIndexSet::whichBetween(TaoVec*V1,TaoVec*V2,TaoVec*V3){
  int i,n,low,high,low2,high2,low3,high3,n_vm=0,info;
  int *vm;
  double *v1,*v2,*vmiddle;

  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoESIIndexSet::getSize"
int TaoESIIndexSet::getSize(int *nn){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#endif
