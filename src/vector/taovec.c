#include "tao_general.h"          /*I "tao_solver.h"  I*/
#include "taovec.h"

#undef __FUNCT__
#define __FUNCT__ "TaoVecDestroy"
/*@C
   TaoVecDestroy - Destroys the TaoVec object.

   Input Parameter:
.  vv - the vector

   Level: beginner
@*/
int TaoVecDestroy( TaoVec* vv){
  TaoFunctionBegin;
  if (vv!=TAO_NULL && vv!=0){ 
    delete vv;
  }
  vv=TAO_NULL;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::CloneVecs"
/*@C
   CloneVecs - Creates an array of pointers to new TaoVec objects. The new 
   objects have the same structure as this one.

   Input Parameter:
.  nn -  number of new vectors

   Output Parameter:
.  tvs -  pointer to array TaoVec pointers

.seealso TaoVec::Clone()

   Level: intermediate
@*/
int TaoVec::CloneVecs(TaoInt nn, TaoVec***tvs){
  int info;
  TaoInt i;
  TaoVec ** ntv;
  TaoFunctionBegin;
  ntv = new TaoVec*[nn];
  for (i=0;i<nn;i++){
    info = this->Clone(&ntv[i]);CHKERRQ(info);
  }
  *tvs=ntv;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::DestroyVecs"
/*@C
   DestroyVecs - Destroys an array TaoVec objects.

   Input Parameter:
+  nn -  number of new vectors
-  tvs -  pointer to array TaoVec pointers

   Level: advanced

.seealso TaoVec::CloneVecs()
@*/
int TaoVec::DestroyVecs(TaoInt nn, TaoVec**ntv){
  int info;
  TaoInt i;
  TaoFunctionBegin;
  for (i=0;i<nn;i++){
    info = TaoVecDestroy( ntv[i] );CHKERRQ(info);
  }
  delete[] ntv;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Clone"
/*@C
   Clone - Creates a new TaoVec object with the same structure as this one.  It does not
   copy the value to the new vector.

   Input:
.  vv - address of a pointer to a TaoVec

   Output Parameter:
.  vv - address of a pointer to new TaoVec object

.seealso TaoVec::CloneVecs(), TaoVec::CopyFrom()

   Level: intermediate
@*/
int TaoVec::Clone( TaoVec* *vv ){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Compatible"
/*@C
   Compatible - Determines whether this vector belongs to the same space as another,
   and operations such as inner product and sum are well defined.

   Input Parameter:
.  vv -  TAO vector to which to the comparison is made

   Output Value:
.  flag - TAO_TRUE if the two vectors are Compatible and TAO_FALSE otherwise.

   Level: advanced

.seealso TaoVec::GetDimension()
@*/
int TaoVec::Compatible(TaoVec* vv, TaoTruth *flag){
  TaoFunctionBegin;
  if (!flag){
    TaoFunctionReturn(1);
  }
  *flag=TAO_FALSE;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::SetToConstant"
/*@C
   SetToConstant - Sets each element of this vector equal to the specified constant.

   Input Parameter:
.  c -  a constant

   Level: intermediate

.seealso TaoVec::Scale()
@*/
int TaoVec::SetToConstant( double c ){
  int info;
  TaoInt nn,i;
  TaoScalar *tptr;

  TaoFunctionBegin;
  info = this->GetArray(&tptr,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ tptr[i]=c; }
  info = this->RestoreArray(&tptr,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::SetToZero"
/*@C
   SetToZero - Sets each element of this vector equal to zero.

   Input Parameters: none

   Level: intermediate

.seealso TaoVec::SetToConstant()
@*/
int TaoVec::SetToZero(){
  int info;
  TaoInt nn,i;
  TaoScalar *tptr;

  TaoFunctionBegin;
  info = this->GetArray(&tptr,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ tptr[i]=0; }
  info = this->RestoreArray(&tptr,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::CopyFrom"
/*@C
   CopyFrom - Copies the contents of one vector into this vector.

   Input Parameter:
.  vv -  A TaoVec from which the contents will be copied.

   Level: intermediate

.seealso TaoVec::Axpy(), TaoVec::ScaleCopyFrom()
@*/
int TaoVec::CopyFrom( TaoVec* vv ){
  int info;
  TaoInt nn1,nn2,i;
  TaoScalar *tptr1,*tptr2;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  if (vv!=this){
    info = vv->GetArray(&tptr2,&nn2);CHKERRQ(info);
    if (nn1!=nn2) {TaoFunctionReturn(1);}
    for (i=0;i<nn1;i++){ tptr1[i]=tptr2[i]; }
    info = vv->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  }
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::ScaleCopyFrom"
/*@C
   ScaleCopyFrom - Copies the contents of one vector into this vector and scales it.

   Input Parameter:
+  a - the scalar 
-  vv -  A TaoVec from which the contents will be copied.

   Level: intermediate

.seealso TaoVec::Axpy(), TaoVec::Aypx()
@*/
int TaoVec::ScaleCopyFrom( double a, TaoVec* vv ){
  int info;
  TaoInt nn1,nn2,i;
  TaoScalar *tptr1,*tptr2;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  if (vv!=this){
    info = vv->GetArray(&tptr2,&nn2);CHKERRQ(info);
    if (nn1!=nn2) {TaoFunctionReturn(1);}
    for (i=0;i<nn1;i++){ tptr1[i]=a*tptr2[i]; }
    info = vv->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  } else {
    this->Scale(a);
  }
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::NormInfinity"
/*@C
   NormInfinity - Computes the infinity norm of this vector.

   Ouput Parameter:
.  vnorm -  the infinity norm of this vector

   Level: intermediate

.seealso TaoVec::Norm1(), TaoVec::Norm2()
@*/
int TaoVec::NormInfinity(double *vnorm){
  int  info;
  TaoInt nn,i;
  TaoScalar dd=0, *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){
    if (v[i]<0 && v[i]<-dd) dd=-v[i];
    else if (v[i]>0 && v[i]>dd) dd=v[i];
  }
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  *vnorm=dd;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Norm1"
/*@C
   Norm1 - Computes the one-norm of this vector.

   Ouput Parameter:
.  vnorm -  the one-norm of this vector

   Level: intermediate

.seealso TaoVec::NormInfinity(), TaoVec::Norm2()
@*/
int TaoVec::Norm1(double *vnorm){
  int info;
  TaoInt nn,i;
  TaoScalar dd=0, *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (v[i]<0) dd-=v[i]; 
    else dd+=v[i];
  }
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  *vnorm = dd;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Norm2"
/*@C
   Norm2 - Compute the two-norm of this vector.

   Ouput Parameter:
.  vnorm -  the two-norm of this vector

   Level: intermediate

.seealso TaoVec::Norm1(), TaoVec::NormInfinity(), TaoVec::Norm2squared()
@*/
int TaoVec::Norm2(double *vnorm){
  int info;
  TaoInt nn,i;
  TaoScalar dd=0, *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) dd+=v[i]*v[i];
  *vnorm=sqrt(dd);
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Norm2squared"
/*@C
   Norm2squared - Computes the square of the two norm of this vector.

   Ouput Parameter:
.  vnorm2 -  the square of the two norm of this vector

   Level: intermediate

.seealso TaoVec::Norm2()
@*/
int TaoVec::Norm2squared(double *vnorm2){
  int info;
  TaoInt nn,i;
  TaoScalar dd=0, *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) dd+=v[i]*v[i];
  *vnorm2=dd;
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Scale"
/*@C
   Scale - Multiplies this vector by a scalar.

   Input Parameter:
.  alpha -  the scalar

   Level: intermediate

.seealso TaoVec::SetToConstant(),  TaoVec::Aypx()
@*/
int TaoVec::Scale( double alpha ){
  int info;
  TaoInt nn,i;
  TaoScalar *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) v[i]*=alpha;
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Axpy"
/*@C
   Axpy - Adds a scalar multiple of a vector to this vector. (this += alpha * xx)

   Input Parameter:
+  alpha -  the scalar 
-  xx - the vector

   Level: intermediate

.seealso TaoVec::CopyFrom(), TaoVec::Aypx() 
@*/
int TaoVec::Axpy( double alpha, TaoVec* xx ){
  int info;
  TaoInt nn1,nn2,i;
  TaoScalar *tptr1,*tptr2;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  info = xx->GetArray(&tptr2,&nn2);CHKERRQ(info);
  if (nn1!=nn2) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++){ tptr1[i]+= alpha * tptr2[i]; }
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  info = xx->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Axpby"
/*@C
   Axpby - Adds a scalar multiple of a vector to a multiple of this vector. (this=alpha*xx + beta*this)

   Input Parameter:
+  alpha -  the scalar of tx
.  xx - the vector
-  beta -  the scalar multiple of this vector

   Level: intermediate

.seealso TaoVec::Axpy(), TaoVec::Aypx() 
@*/
int TaoVec::Axpby( double alpha, TaoVec* xx, double beta ){
  int info;
  TaoInt nn1,nn2,i;
  TaoScalar *tptr1,*tptr2;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  info = xx->GetArray(&tptr2,&nn2);CHKERRQ(info);
  if (nn1!=nn2) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++){ tptr1[i] = beta * tptr1[i] + alpha * tptr2[i]; }
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  info = xx->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Aypx"
/*@C
   Aypx - Adds a vector to a scalar multiple of this vector. (this=alpha*this+xx)

   Input Parameter:
+  alpha -  the scalar 
-  xx - the vector

   Level: intermediate

.seealso TaoVec::Scale(), TaoVec::Axpy()
@*/
int TaoVec::Aypx( double alpha, TaoVec* xx ){
  int info;
  TaoInt nn1,nn2,i;
  TaoScalar *tptr1,*tptr2;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  info = xx->GetArray(&tptr2,&nn2);CHKERRQ(info);
  if (nn1!=nn2) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++){ tptr1[i] = alpha * tptr1[i] + tptr2[i]; }
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  info = xx->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  TaoFunctionReturn(0);
}
 
#undef __FUNCT__
#define __FUNCT__ "TaoVec::AddConstant"
/*@C
   AddConstant - Adds a constant to each element of this vector.

   Input Parameter:
.  alpha -  the scalar 

   Level: intermediate

.seealso TaoVec::SetToConstant(), TaoVec::Axpy() 
@*/
int TaoVec::AddConstant( double alpha ){
  int info;
  TaoInt nn,i;
  TaoScalar *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) v[i]+=alpha;
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Dot"
/*@C
   Dot - Computes the inner product of this vector with another vector.

   Input Parameter:
.  vv -  another TaoVec object

   Output Parameter:
.  vDotv -  the inner product of the two vectors

   Level: intermediate

.seealso TaoVec::Norm() 
@*/
int TaoVec::Dot( TaoVec* vv, double *vDotv ){
  int info;
  TaoInt nn1,nn2,i;
  TaoScalar dd=0,*tptr1,*tptr2;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  info = vv->GetArray(&tptr2,&nn2);CHKERRQ(info);
  if (nn1!=nn2) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++) dd+=tptr1[i]*tptr2[i];
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  info = vv->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  *vDotv=dd;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Negate"
/*@C
   Negate - Multiplies the elements of this vector by negative one.

   Input Parameters: none

   Level: intermediate

.seealso TaoVec::Scale()
@*/
int TaoVec::Negate(){ 
  int info;
  TaoInt nn,i;
  TaoScalar *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++) v[i]=-v[i];
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Reciprocal"
/*@C
   Reciprocal - Sets each element of this vector to its Reciprocal.

   Input Parameters: none

   Level: intermediate

.seealso TaoVec::PointwiseDivide()
@*/
int TaoVec::Reciprocal(){ 
  int info;
  TaoInt nn,i;
  TaoScalar *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (v[i]!=0)  v[i]= 1.0/v[i];
    else v[i]=TAO_INFINITY;
  }
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Sqrt"
/*@C
   Sqrt - Sets each element of this vector to its square root.

   Input Parameters: none

   Level: intermediate

@*/
int TaoVec::Sqrt(){ 
  int info;
  TaoInt i;
  TaoInt nn;
  TaoScalar *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (v[i] >= 0) {
      v[i] = sqrt(v[i]);
    }
    else {
      v[i] = TAO_INFINITY;
    }
  }
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Pow"
/*@C
   Pow - Raises each element of this vector to a power.

   Input Parameters: none

   Level: intermediate

@*/
int TaoVec::Pow(double p){ 
  int info;
  TaoInt nn,i;
  TaoScalar *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (v[i] >= 0) {
      v[i] = pow(v[i], p);
    }
    else {
      v[i] = TAO_INFINITY;
    }
  }
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::GetDimension"
/*@C
   GetDimension - Gets the dimension of the vector space where this vector belongs.

   Output Parameter:
.  n - the dimension of the vector space

   Level: intermediate

.seealso TaoVec::Compatible()
@*/
int TaoVec::GetDimension(TaoInt *n){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::PointwiseMultiply"
/*@C
   PointwiseMultiply - Computes the componentwise multiplication of two vectors 
   and stores the result in this vector.

   Input Parameters:
.   vv, ww - the two vectors

   Level: intermediate

.seealso TaoVec::PointwiseDivide()
@*/
int TaoVec::PointwiseMultiply( TaoVec* vv, TaoVec* ww ){
  int info;
  TaoInt nn1,nn2,nn3,i;
  TaoScalar *tptr1,*tptr2,*tptr3;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  if (this!=vv){ 
    info = vv->GetArray(&tptr2,&nn2);CHKERRQ(info);
  } else {
    tptr2=tptr1; nn2=nn1;
  }
  if (this!=ww && vv!=ww){ 
    info = ww->GetArray(&tptr3,&nn3);CHKERRQ(info);
  } else if (vv==ww){ 
    tptr3=tptr2; nn3=nn2;
  } else {
    tptr3=tptr1; nn3=nn1;
  }
  if (nn1!=nn2 || nn2!=nn3) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++) tptr1[i]=tptr2[i] * tptr3[i];
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  if (this!=vv){ 
    info = vv->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  }
  if (this!=ww && vv!=ww){ 
    info = ww->RestoreArray(&tptr3,&nn3);CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::PointwiseDivide"
/*@C
   PointwiseDivide - Computes the componentwise division of two vectors 
   and stores the result in this vector.

   Input Parameters:
+  vv - the vector of numerators
-  ww - the vector of denominators

   Level: intermediate

.seealso TaoVec::PointwiseMultiply()
@*/
int TaoVec::PointwiseDivide( TaoVec* vv , TaoVec* ww){
  int info;
  TaoInt nn1,nn2,nn3,i;
  TaoScalar *tptr1,*tptr2,*tptr3;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  if (this!=vv){ 
    info = vv->GetArray(&tptr2,&nn2);CHKERRQ(info);
  } else {
    tptr2=tptr1; nn2=nn1;
  }
  if (this!=ww && vv!=ww){ 
    info = ww->GetArray(&tptr3,&nn3);CHKERRQ(info);
  } else if (vv==ww){
    tptr3=tptr2; nn3=nn2;
  } else {
    tptr3=tptr1; nn3=nn1;
  }
  if (nn1!=nn2 || nn2!=nn3) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++) tptr1[i]=tptr2[i] / tptr3[i];
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  if (this!=vv){ 
    info = vv->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  }
  if (this!=ww && vv!=ww){ 
    info = ww->RestoreArray(&tptr3,&nn3);CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Median"
/*@C
   Median - Computes the componentwise median of three vectors 
   and stores the result in this vector.

   Input Parameters:
.   vv, ww, xx - the three vectors

   Level: intermediate

@*/
int TaoVec::Median( TaoVec* vv, TaoVec* ww, TaoVec* xx){
  int info;
  TaoInt nn1,nn2,nn3,nn4,i;
  TaoScalar *tptr1,*tptr2,*tptr3,*tptr4;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  if (this!=vv){
    info = vv->GetArray(&tptr2,&nn2);CHKERRQ(info);
  } else {
    tptr2=tptr1; nn2=nn1;
  }
  if (this!=ww){
    info = ww->GetArray(&tptr3,&nn3);CHKERRQ(info);
  } else {
    tptr3=tptr1; nn3=nn1;
  }
  if (this!=xx){
    info = xx->GetArray(&tptr4,&nn4);CHKERRQ(info);
  } else {
    tptr4=tptr1; nn4=nn1;
  }

  if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4) {TaoFunctionReturn(1);}
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
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  if (this!=vv){
    info = vv->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  }
  if (this!=ww){
    info = ww->RestoreArray(&tptr3,&nn3);CHKERRQ(info);
  }
  if (this!=xx){
    info = xx->RestoreArray(&tptr4,&nn4);CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::PointwiseMinimum"
/*@C
   PointwiseMinimum - Computes the componentwise minimum of two vectors 
   and stores the result in this vector.

   Input Parameters:
.   vv, ww - the two vectors

   Level: intermediate

.seealso TaoVec::PointwiseMaximum()
@*/
int TaoVec::PointwiseMinimum( TaoVec* vv, TaoVec* ww){
  int info;
  TaoInt nn1,nn2,nn3,i;
  TaoScalar *tptr1,*tptr2,*tptr3;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  if (vv!=this){
    info = vv->GetArray(&tptr2,&nn2);CHKERRQ(info);
  } else {
    tptr2=tptr1;  nn2=nn1;
  }
  if (ww!=this){
    info = ww->GetArray(&tptr3,&nn3);CHKERRQ(info);
  } else {
    tptr3=tptr1;  nn3=nn1;
  }

  if (nn1!=nn2 || nn2!=nn3) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++) tptr1[i] = TaoMin( tptr2[i] , tptr3[i]);
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  if (vv!=this){
    info = vv->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  }
  if (ww!=this){
    info = ww->RestoreArray(&tptr3,&nn3);CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::PointwiseMaximum"
/*@C
   PointwiseMaximum - Computes the componentwise minimum of two vectors 
   and stores the result in this vector.

   Input Parameters:
.  vv, ww - the two vectors

   Level: intermediate

.seealso TaoVec::PointwiseMinimum()
@*/
int TaoVec::PointwiseMaximum( TaoVec* vv, TaoVec* ww){
  int info;
  TaoInt nn1,nn2,nn3,i;
  TaoScalar *tptr1,*tptr2,*tptr3;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  if (vv!=this){
    info = vv->GetArray(&tptr2,&nn2);CHKERRQ(info);
  } else {
    tptr2=tptr1;  nn2=nn1;
  }
  if (ww!=this){
    info = ww->GetArray(&tptr3,&nn3);CHKERRQ(info);
  } else {
    tptr3=tptr1;  nn3=nn1;
  }
  if (nn1!=nn2 || nn2!=nn3) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++) tptr1[i] = TaoMax( tptr2[i] , tptr3[i]);
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  if (vv!=this){
    info = vv->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  }
  if (ww!=this){
    info = ww->RestoreArray(&tptr3,&nn3);CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Waxpby"
/*@C
   Waxpby - Sums two scaled vectors and stores the result in this vector. (this=alpha*xx+beta*yy)

   Input Parameters:
+  a - the multiple of the first vector
.  xx - the first vector
.  b - the multiple of the second vector
-  yy - the second vector

   Level: intermediate

.seealso TaoVec::Axpy(), TaoVec::Axpby(), TaoVec::Aypx();
@*/
int TaoVec::Waxpby  ( double a, TaoVec* xx, double b, TaoVec* yy){
  int info;
  TaoInt nn1,nn2,nn3,i;
  TaoScalar *tptr1,*tptr2,*tptr3;

  TaoFunctionBegin;
  info = this->GetArray(&tptr1,&nn1);CHKERRQ(info);
  info = xx->GetArray(&tptr2,&nn2);CHKERRQ(info);
  info = yy->GetArray(&tptr3,&nn3);CHKERRQ(info);
  if (nn1!=nn2 || nn2!=nn3) {TaoFunctionReturn(1);}
  for (i=0;i<nn1;i++){ tptr1[i] = a * tptr2[i] + b * tptr3[i]; }
  info = this->RestoreArray(&tptr1,&nn1);CHKERRQ(info);
  info = xx->RestoreArray(&tptr2,&nn2);CHKERRQ(info);
  info = yy->RestoreArray(&tptr3,&nn3);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::AbsoluteValue"
/*@C
   AbsoluteValue - Sets each element of this vector equal to its magnitude.

   Input Parameters: none

   Level: intermediate
@*/
int TaoVec::AbsoluteValue(){
  int info;
  TaoInt nn,i;
  TaoScalar *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    v[i]= TaoAbsScalar(v[i]);
  }
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::MinElement"
/*@C
   MinElement - Finds the smallest element of this vector.

   Output Parameter:
.  val - the smallest value in the vector

   Level: intermediate
@*/
int TaoVec::MinElement(double *val){
  int info;
  TaoInt nn,i;
  TaoScalar dd=TAO_INFINITY,*v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){ 
    if (v[i]<dd) dd=v[i];
  }
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  *val = dd;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Sign"
/*@C
   Sign - Replace each element of this vector with -1, 0, or 1, depending on its sign.

   Level: intermediate
@*/
int TaoVec::Sign(){
  int info;
  TaoInt nn,i;
  TaoScalar *v;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){
    if (v[i]<0){ 
      v[i]=-1;
    } else if (v[i]>0){
      v[i]=1.0;
    } else {
      v[i]=0.0;
    }
  }
  info = this->RestoreArray(&v,&nn);CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::SetReducedVec"
/*@C
   SetReducedVec - Sets the reduced space of the vector that this
   vector should represent.  The index set describes which
   elements of the vector should be used.  This routine also
   copies  the appropriate elements from the full space vector
   into the reduced space vector.

   Input Parameters:
+  vv -  a vector
-  ss -  an index set

   Level: advanced

.seealso TaoVec::CreateReducedVec(),  TaoVec::ReducedCopyFromFull(), TaoVec::ReducedXPY(), TaoSelectSubset()
@*/
int TaoVec::SetReducedVec(TaoVec* vv,TaoIndexSet* ss){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::ReducedCopyFromFull"
/*@C
  ReducedCopyFromFull - Copies the appropriate elements of the vector into
  this reduced vector.
  
  Input Parameters:
+  ss -  an index set
-  vv -  a full vector

   Level: advanced

.seealso TaoVec::CreateReducedVec(), TaoVec::ReducedXPY()
@*/
int TaoVec::ReducedCopyFromFull(TaoVec* vv, TaoIndexSet* ss){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::ReducedXPY"
/*@C
  ReducedXPY - Adds a reduced vector to the appropriate elements of this vector.
  
  Input Parameters:
+  vv -  the reduced vector
-  ss -  the index set identifying which elements of this vector should be supplemented

   Level: advanced

.seealso TaoVec::CreateReducedVec(), TaoVec::ReducedCopyFromFull()
@*/
int TaoVec::ReducedXPY(TaoVec* vv, TaoIndexSet* ss){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::StepMax"
/*@C
  StepMax - Calculates the largest multiple of a vector that can be added to
  this vector while keeping each element of this vector nonnegative.
  
  Input Parameters:
. vv - the step direction

  Input Parameters:
. step1 - the maximum stepsize

  Note:
  This vector should contain all positive elements.

  Note:
  If there is no maximum steplength, the output scalar may be set
  to TAO_INFINITY.

  Level: advanced
@*/
int TaoVec::StepMax(TaoVec* vv,double *step1){
  int info;
  TaoInt nn1,nn2,i;
  TaoScalar *xx,*dx;
  double stepmax1=TAO_INFINITY;

  TaoFunctionBegin;
  info = this->GetArray(&xx,&nn1);CHKERRQ(info);
  if (this!=vv){
    info = vv->GetArray(&dx,&nn2);CHKERRQ(info);
    if (nn1!=nn2) {TaoFunctionReturn(1);}
    for (i=0;i<nn1;i++){
      if (xx[i] < 0){
	TaoFunctionReturn(1);
      } else if (dx[i]<0){ 
	stepmax1=TaoMin(stepmax1,-xx[i]/dx[i]);
      }
    }
    info = vv->RestoreArray(&dx,&nn2);CHKERRQ(info);
    *step1=stepmax1;
  } else {
    *step1=1.0;
  }
  info = this->RestoreArray(&xx,&nn1);CHKERRQ(info);
  TaoFunctionReturn(0);
  }

/*@C
  StepBoundInfo - Calculates the largest multiple of a vector that can be added to
  this vector while keeping each element of this vector nonnegative.
  
  Input Parameters:
+ txl - the lower bounds on this vector
. txu - the upper bounds on this vector
- tdx - the step direction for this vector

  Output Parameters:
+ boundmin - the step to closest bound i.e min(a1, ..., an);
. wolfemin - the step to closest bound not equal i.e min(b1, ..., bn);
- boundmax - the step to farthest bound   i.e. max(c1, ..., cn);

  Where:
  if tdx[i] > 0;  ai = (txu[i] - this[i])/tdx[i] ; bi=ai, ci=ai;
  if tdx[i] < 0;  ai = (txl[i] - this[i])/tdx[i] ; bi=ai, ci=ai
  if tdx[i] == 0 && txl[i] < x[i] < txu[i] ; ai=TAO_INFINITY, bi=ai, ci=ai;
  if tdx[i] == 0 && (txl[i] == x[i] || txu[i] == x[i]) ; ai= 0, bi=TAO_INFINITY, ci=0;
 
  Note:
  If there is no maximum steplength, the output scalar may be set
  to TAO_INFINITY.

  Level: advanced
@*/
int TaoVec::StepBoundInfo(TaoVec* XL ,TaoVec* XU, TaoVec*S, double *bmin1,double *bmin2, double *bmax){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::View"
/*@C
  View - Views the contents of the vector.
  
  Input Parameters: none
  
  Level: intermediate
@*/
int TaoVec::View(){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::BoundGradientProjection"
/*@C
  BoundGradientProjection - Projects this vector according to this definition.
  If XX[i]==XXL[i], then this[i] = min(G[i],0);   
  If XX[i]==XXU[i], then this[i] = max(G[i],0);
  else this[i] = G[i];

  Input Parameters:
. GG,XX,XXL,XXU - the vectors.

  Level: advanced
@*/
int TaoVec::BoundGradientProjection(TaoVec* gg,TaoVec* xxll,TaoVec* xx, TaoVec* xxuu){
  int info;
  TaoInt nn1,nn2,nn3,nn4,nn5,i;
  TaoScalar *xptr,*xlptr,*xuptr,*gptr,*gpptr;

  TaoFunctionBegin;
  info = this->GetArray(&gpptr,&nn1);CHKERRQ(info);
  if (this != gg){
    info = gg->GetArray(&gptr,&nn2);CHKERRQ(info);
  } else {
    gptr=gpptr; nn2=nn1;
  }
  info = xxll->GetArray(&xlptr,&nn3);CHKERRQ(info);
  info = xx->GetArray(&xptr,&nn4);CHKERRQ(info);
  info = xxuu->GetArray(&xuptr,&nn5);CHKERRQ(info);
  if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4 || nn4!=nn5) {TaoFunctionReturn(1);}

  for (i=0; i<nn1; i++){

    gpptr[i] = gptr[i];
    if (gpptr[i]>0 && xptr[i]<=xlptr[i]){
      gpptr[i] = 0;
    } else if (gpptr[i]<0 && xptr[i]>=xuptr[i]){
      gpptr[i] = 0;
    }
  }
  info = this->RestoreArray(&gpptr,&nn1);CHKERRQ(info);
  if (this!=gg){
    info = gg->RestoreArray(&gptr,&nn2);CHKERRQ(info);
  }
  info = xxll->RestoreArray(&xlptr,&nn3);CHKERRQ(info);
  info = xx->RestoreArray(&xptr,&nn4);CHKERRQ(info);
  info = xxuu->RestoreArray(&xuptr,&nn5);CHKERRQ(info);

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::GetArray"
/*@C
   GetArray - Sets a pointer to the first element in the vector array.

   Output Parameters:
+  dptr -  pointer an the array of numbers
-  n -  the length of the array

   Note:
   This operation may not be defined the same for all vector types.

   Level: intermediate

.seealso TaoVec::RestoreArray()
@*/
int TaoVec::GetArray(TaoScalar **dptr, TaoInt *n){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::RestoreArray"
/*@C
   RestoreArray - Returns a pointer to the first element in the vector array.

   Input Parameters:
+  dptr -  pointer an the array of numbers obtained from TaoVec::GetArray
-  n -  the length of the array

   Note:
   This routine is not used within TAO solvers.  Rather, it is used to
   implement other routines.

   Level: intermediate

.seealso TaoVec::GetArray()
@*/
int TaoVec::RestoreArray(TaoScalar **dptr, TaoInt *n){
  TaoFunctionBegin;
  *dptr=NULL;
  *n=0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::CreateIndexSet"
/*@C
   CreateIndexSet - Creates an index set that may be used to describe sets of
   elements of this vector.

   Output Parameters:
.  ss -  a new index set

   Level: advanced
@*/
int TaoVec::CreateIndexSet(TaoIndexSet **ss){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

inline static double fischer(double a, double b) 
{
   // Method suggested by Bob Vanderbei
   if (a + b <= 0) {
     return sqrt(a*a + b*b) - (a + b);
   }
   return -2.0*a*b / (sqrt(a*a + b*b) + (a + b));
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::Fischer"
/*@C
   Fischer - Evaluates the Fischer-Burmeister function for complementarity 
   problems.
   
   Input Parameters:
+  xx - current point
.  ff - function evaluated at x
.  ll - lower bounds 
-  uu - upper bounds

   Notes: 
   The Fischer-Burmeister function is defined as
$        phi(a,b) := sqrt(a*a + b*b) - a - b
   and is used reformulate a complementarity problem as a semismooth
   system of equations.

   The result of this function is done by cases:
+  l[i] == -infinity, u[i] == infinity  -- fb[i] = -f[i]
.  l[i] == -infinity, u[i] finite       -- fb[i] = phi(u[i]-x[i], -f[i])
.  l[i] finite,       u[i] == infinity  -- fb[i] = phi(x[i]-l[i],  f[i])
.  l[i] finite < u[i] finite -- fb[i] = phi(x[i]-l[i], phi(u[i]-x[i], -f[u]))
-  otherwise l[i] == u[i] -- fb[i] = l[i] - x[i]

   Level: advanced

.seealso TaoMat::D_Fischer()
@*/
int TaoVec::Fischer(TaoVec* xx, TaoVec* ff, TaoVec* ll, TaoVec* uu)
{
  int info;
  TaoInt nn1,nn2,nn3,nn4,nn5,i;
  TaoScalar *v,*x,*f,*l,*u;

  TaoFunctionBegin;
  info = this->GetArray(&v,&nn1); CHKERRQ(info);
  info = xx->GetArray(&x,&nn2); CHKERRQ(info);
  info = ff->GetArray(&f,&nn3); CHKERRQ(info);
  info = ll->GetArray(&l,&nn4); CHKERRQ(info);
  info = uu->GetArray(&u,&nn5); CHKERRQ(info);
  if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4 || nn4!=nn5) { 
    TaoFunctionReturn(1); 
  }

  for (i=0;i<nn1;i++) {

    if ((l[i] <= -TAO_INFINITY) && (u[i] >= TAO_INFINITY)) {
      v[i] = -f[i];
    } 
    else if (l[i] <= -TAO_INFINITY) {
      v[i] = -fischer(u[i] - x[i], -f[i]);
    } 
    else if (u[i] >=  TAO_INFINITY) {
      v[i] =  fischer(x[i] - l[i],  f[i]);
    } 
    else if (l[i] == u[i]) {
      v[i] = l[i] - x[i];
    } 
    else {
      v[i] =  fischer(u[i] - x[i], -f[i]);
      v[i] =  fischer(x[i] - l[i],  v[i]);
    }

  }

  info = this->RestoreArray(&v,&nn1);CHKERRQ(info);
  info = xx->RestoreArray(&x,&nn2);CHKERRQ(info);
  info = ff->RestoreArray(&f,&nn3);CHKERRQ(info);
  info = ll->RestoreArray(&l,&nn4);CHKERRQ(info);
  info = uu->RestoreArray(&u,&nn5);CHKERRQ(info);

  TaoFunctionReturn(0);
}

inline static double sfischer(double a, double b, double c)
{
   // Method suggested by Bob Vanderbei
   if (a + b <= 0) {
     return sqrt(a*a + b*b + 2.0*c*c) - (a + b);
   }
   return 2.0*(c*c - a*b) / (sqrt(a*a + b*b + 2.0*c*c) + (a + b));
}

#undef __FUNCT__
#define __FUNCT__ "TaoVec::SFischer"
/*@C
   SFischer - Evaluates the Smoothed Fischer-Burmeister function for
   complementarity problems.

   Input Parameters:
+  xx - current point
.  ff - function evaluated at x
.  ll - lower bounds
.  uu - upper bounds
-  mu - smoothing parameter

   Notes:
   The Smoothed Fischer-Burmeister function is defined as
$        phi(a,b) := sqrt(a*a + b*b + 2*mu*mu) - a - b
   and is used reformulate a complementarity problem as a semismooth
   system of equations.

   The result of this function is done by cases:
+  l[i] == -infinity, u[i] == infinity  -- fb[i] = -f[i] - 2*mu*x[i]
.  l[i] == -infinity, u[i] finite       -- fb[i] = phi(u[i]-x[i], -f[i], mu)
.  l[i] finite,       u[i] == infinity  -- fb[i] = phi(x[i]-l[i],  f[i], mu)
.  l[i] finite < u[i] finite -- fb[i] = phi(x[i]-l[i], phi(u[i]-x[i], -f[u], mu), mu)
-  otherwise l[i] == u[i] -- fb[i] = l[i] - x[i]

   Level: advanced

.seealso TaoMat::SD_Fischer()
@*/
int TaoVec::SFischer(TaoVec* xx, TaoVec* ff, TaoVec* ll, TaoVec* uu, double mu)
{

  int info;
  TaoInt nn1, nn2, nn3, nn4, nn5,i;
  TaoScalar *v, *x, *f, *l, *u;

  TaoFunctionBegin;

  if ((mu >= -TAO_EPSILON) && (mu <= TAO_EPSILON)) {
    Fischer(xx, ff, ll, uu);
  }
  else {
    info = this->GetArray(&v, &nn1); CHKERRQ(info);
    info = xx->GetArray(&x, &nn2); CHKERRQ(info);
    info = ff->GetArray(&f, &nn3); CHKERRQ(info);
    info = ll->GetArray(&l, &nn4); CHKERRQ(info);
    info = uu->GetArray(&u, &nn5); CHKERRQ(info);

    if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4 || nn4!=nn5) {
      TaoFunctionReturn(1);
    }

    for (i = 0; i < nn1; ++i) {
      if ((l[i] <= -TAO_INFINITY) && (u[i] >= TAO_INFINITY)) {
        v[i] = -f[i] - mu*x[i];
      } 
      else if (l[i] <= -TAO_INFINITY) {
        v[i] = -sfischer(u[i] - x[i], -f[i], mu);
      } 
      else if (u[i] >=  TAO_INFINITY) {
        v[i] =  sfischer(x[i] - l[i],  f[i], mu);
      } 
      else if (l[i] == u[i]) {
        v[i] = l[i] - x[i];
      } 
      else {
        v[i] =  sfischer(u[i] - x[i], -f[i], mu);
        v[i] =  sfischer(x[i] - l[i],  v[i], mu);
      }
    }

    info = this->RestoreArray(&v, &nn1); CHKERRQ(info);
    info = xx->RestoreArray(&x, &nn2); CHKERRQ(info);
    info = ff->RestoreArray(&f, &nn3); CHKERRQ(info);
    info = ll->RestoreArray(&l, &nn4); CHKERRQ(info);
    info = uu->RestoreArray(&u, &nn5); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

