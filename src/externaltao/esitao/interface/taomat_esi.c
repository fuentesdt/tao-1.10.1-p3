#include "esi/ESI.h"
#include "taomat_esi.h"
#include "tao_general.h"


#undef __FUNCT__
#define __FUNCT__ "TaoWrapESIMat"
/*@C
   TaoWrapESIMat - Creates a new TaoMat object using a PETSc matrix.

   Input Parameter:
+  M -  a PETSc matrix
.  flg - set to TAO_TRUE if the PETSc matrtix M should be destroyed when this TaoMat is destroyed, and set to TAO_FALSE otherwise.
-  MM - new TaoMat

   Level: beginner

.seealso TaoMatGetPetscMat(), TaoMatDestroy()
@*/
int TaoWrapESIMat(esi::Operator<double, int> *M, TaoMatESI ** MM){
  TaoFunctionBegin;
  *MM = new  TaoMatESI(M);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatGetESIMat"
/*@C
   TaoMatGetESIMat - If the TaoMat is of the TaoMatPetsc type, this routine
   gets the underlying PETSc matrix.

   Input Parameter:
.  MM - the TaoMatESI 

   Output Parameter:
.  M -  the PETSc mat

   Note:
   The function TaoMatESI::GetMat() will also return the Mat M.

   Level: beginner
@*/
int TaoMatGetESIMat( TaoMat *MM, esi::Operator<double, int> ** M){
  TaoFunctionBegin;
  *M=((TaoMatESI *)MM)->esiobj;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::TaoMatESI"
/* @C
   TaoMatESI - Creates a new TaoMat object using a PETSc matrix.

   Input Parameter:
+  MM -  a PETSc matrix
-  flg - set to TAO_TRUE if VV should be destroyed when it is replaced or when this TaoMat is deleted.  If set to TAO_FALSE, the destruction of VV must
be done elsewhere.

   Level: beginner


   Notes:
   The method TaoMatESI::GetMat() returns a pointer to the matrix MM

   The method TaoMatESI::replaceMat(Mat M) replaces the matrix MM with M.
@ */
TaoMatESI::TaoMatESI( esi::Operator<double, int>* MM)
  :TaoMat(){
  esiobj=MM;
  //  esiobj2=MM;
  MM->addReference();
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::Compatible"
int TaoMatESI::Compatible(TaoVec *xx, TaoVec *yy, TaoTruth *flag){
  *flag = TAO_TRUE;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::GetDimensions"
int TaoMatESI::GetDimensions( int*m, int*n ){
  esi::ErrorCode info;
  esi::MatrixData<int>* esimatwrite=NULL;
  TaoFunctionBegin;
  esimatwrite = dynamic_cast<esi::MatrixData<int>*>(esiobj);
  if (esimatwrite==NULL){
    SETERRQ(56,"Operation not defined");
  } else {
    info=esimatwrite->getGlobalSizes(*m,*n); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::Multiply"
int TaoMatESI::Multiply(TaoVec*tv,TaoVec*tw){
  esi::Vector<double,int> *vv=(esi::Vector<double,int>*)(tv->VecObject);
  esi::Vector<double,int> *ww=(esi::Vector<double,int>*)(tw->VecObject);
  esi::ErrorCode info;
  TaoFunctionBegin;
  info=esiobj->apply(*vv,*ww); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::MultiplyTranspose"
int TaoMatESI::MultiplyTranspose(TaoVec*tv,TaoVec*tw){
  esi::Vector<double,int> *vv=(esi::Vector<double,int>*)(tv->VecObject);
  esi::Vector<double,int> *ww=(esi::Vector<double,int>*)(tw->VecObject);
  esi::OperatorTranspose<double,int>* esioptranspose;
  esi::ErrorCode info;
  TaoFunctionBegin;
  esioptranspose = dynamic_cast<esi::OperatorTranspose<double,int>*>(esiobj);
  if (esioptranspose == NULL){
    SETERRQ(56,"Operation not defined");
  } else {
    info=esioptranspose->applyTranspose(*vv,*ww); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}


int ESIaddDiagonal(esi::Object *,  esi::Vector<double,int> *);
#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::AddDiagonal"
int TaoMatESI::AddDiagonal(TaoVec*tv){
  esi::Vector<double,int> *x=(esi::Vector<double,int>*)(tv->VecObject);
  TaoVecESI * xx = dynamic_cast<TaoVecESI*>(tv);
  x= xx->esivec;

  TaoFunctionBegin;
  esi::ErrorCode info=ESIaddDiagonal(esiobj,x); CHKERRQ(info);
  TaoFunctionReturn(0);
}



int ESIsetDiagonal(esi::Object *,  esi::Vector<double,int> *);
#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::SetDiagonal"
int TaoMatESI::SetDiagonal(TaoVec*tv){
  esi::Vector<double,int> *x=(esi::Vector<double,int>*)(tv->VecObject);
  TaoVecESI * xx = dynamic_cast<TaoVecESI*>(tv);
  x= xx->esivec;

  TaoFunctionBegin;
  esi::ErrorCode info=ESIsetDiagonal(esiobj,x); CHKERRQ(info);
  TaoFunctionReturn(0);
}

int ESIgetDiagonal(esi::Object *,  esi::Vector<double,int> *);
#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::GetDiagonal"
int TaoMatESI::GetDiagonal(TaoVec*tv){
  esi::Vector<double,int> *x=(esi::Vector<double,int>*)(tv->VecObject);
  TaoVecESI * xx = dynamic_cast<TaoVecESI*>(tv);
  x= xx->esivec;

  TaoFunctionBegin;
  esi::ErrorCode info=ESIgetDiagonal(esiobj,x); CHKERRQ(info);
  TaoFunctionReturn(0);
}


int ESIshiftDiagonal( esi::MatrixRowWriteAccess<double,int> *, double);
#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::ShiftDiagonal"
int TaoMatESI::ShiftDiagonal(TaoScalar c){
  esi::ErrorCode info;
  TaoFunctionBegin;
  esi::MatrixRowWriteAccess<double,int>* esimatwrite;
  esimatwrite = dynamic_cast<esi::MatrixRowWriteAccess<double,int>*>(esiobj);
  if (esimatwrite==NULL){
    SETERRQ(56,"Operation not defined");
  } else {
    info = ESIshiftDiagonal(esimatwrite,c); CHKERRQ(info); 
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::norm1"
int TaoMatESI::Norm1(double *nm){
  TaoFunctionBegin;
  *nm=1.0;
  TaoFunctionReturn(0);
}

esi::ErrorCode ESIRowScale( esi::Object *, esi::Vector<double,int>*);
#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::ColScale"
int TaoMatESI::RowScale(TaoVec* tv){
  esi::Vector<double,int> *scale=(esi::Vector<double,int>*)(tv->VecObject);
  TaoFunctionBegin;
  esi::ErrorCode info=ESIRowScale(this->esiobj,scale); CHKERRQ(info);
  TaoFunctionReturn(0);
};

int TaoPrint_ESI_RowMatrix(esi::Object * esiobj);
#undef __FUNCT__
#define __FUNCT__ "TaoMatESI::View"
int TaoMatESI::View(){

  TaoFunctionBegin;
  int info=TaoPrint_ESI_RowMatrix(this->esiobj); CHKERRQ(info);
  TaoFunctionReturn(0);
}

