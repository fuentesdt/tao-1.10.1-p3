/**************************************************************
File: taovec_ga.c

TAO/GLobal Array Project

Author: Limin Zhang, Ph.D.
        Mathematics Department
        Columbia Basin College
        Pasco, WA 99301
        Limin.Zhang@cbc2.org

Mentor: Jarek Naplocha, Ph.D.
        Environmental Molecular Science Laboratory
        Pacific Northwest National Laboratory
        Richland, WA 99352

Purpose:
      to design and implement some interfaces between TAO and
      global arrays.

Date: 4/22/2002


Revise history: 

7/15/02
	 replace pv with this->pv 
	 introduce GAVec as int 

8/8/02
	 clean Petsc function calls and marcos.

**************************************************************/

#include "math.h"


#include "tao_general.h"
#include "taovec_ga.h" /*I "taovec_ga.h" I*/


#undef __FUNCT__
#define __FUNCT__ "TaoWrapGaVec"
/*@C 
  TaoWrapGaVec - Creates a new TaoVec object using an existing
  GlobalArray vector.

   Input Parameter:
+  V -  a GlobalArray vector
-  TV -  the address of a pointer to a TaoGaVec

   Output Parameter:
.  TV - pointer to a new TaoGaVec

   Level: advanced

.seealso TaoVecGetGaVec(), TaoWrapPetscVec()
@*/
int TaoWrapGaVec (GAVec V, TaoVecGa ** TV)
{
  TaoFunctionBegin;
  *TV = new TaoVecGa (V);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGetGaVec"
/*@C
   TaoVecGetGaVec - If the TaoVec is of the TaoVecGa type, this routine gets the underlying Ga Vec.

   Input Parameter:
.  TV - the TaoVecGa 

   Output Parameter:
.  V -  the Ga vector


   Level: advanced

.seealso TaoWrapGaVec()

@*/

int TaoVecGetGaVec (TaoVec * TV, GAVec * V)
{
  TaoFunctionBegin;
  if (TV)
    {
      TaoVecGa *tmp = dynamic_cast<TaoVecGa *>(TV);
      if (!tmp) {
	SETERRQ(1,"Could not cast from TaoVec* to TaoVecGa*.");
      }
      *V = tmp->GetVec ();
    }
  else
    {
      *V = 0;
    }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::TaoVecGa"
TaoVecGa::TaoVecGa (GAVec VV)
{
  this->pv = VV;
  //this->VecObject = (void *) (&VV);
  return;
}

#undef __FUNCT__ 
#define __FUNCT__ "TaoVecGa::Clone"
int
TaoVecGa::Clone (TaoVec ** ntv)
{
  TaoVecGa *nptv;
  GAVec npv;

  TaoFunctionBegin;
  npv = GA_Duplicate (this->pv, (char*)"Cloned");	//pv is the ga handle of the current object. 
  if (!npv)
    GA_Error ((char*)"TaoVecGa::Clone:duplicate failed: ", npv);
  TaoWrapGaVec (npv, &nptv);
  *ntv = nptv;
  npv = 0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Compatible"
int TaoVecGa::Compatible (TaoVec * tv, TaoTruth *flag)
{
  int info;
  int  pv2 = 0;

  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv2); CHKERRQ(info);
  info = GA_Compare_distr (this->pv, pv2);

  if (info == 0)
    *flag =  TAO_TRUE;
  else
    *flag = TAO_FALSE;

  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::SetToConstant"
int
TaoVecGa::SetToConstant (TaoScalar c)
{
  TaoFunctionBegin;
  GA_Fill (this->pv, &c);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::SetToZero"
int
TaoVecGa::SetToZero ()
{
  TaoFunctionBegin;
  GA_Zero (this->pv);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::CopyFrom"
int
TaoVecGa::CopyFrom (TaoVec * tv)
{
  int pv2 = 0;
  TaoFunctionBegin;
  TaoVecGetGaVec (tv, &pv2);
  GA_Copy (pv2, this->pv);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::ScaleCopyFrom"
int
TaoVecGa::ScaleCopyFrom (TaoScalar a, TaoVec * tv)
{
  int info, pv2 = 0;
  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv2);
  GA_Copy (pv2, this->pv);
  GA_Scale (this->pv, &a);
  TaoFunctionReturn(0);
}



#undef __FUNCT__ 
#define __FUNCT__ "TaoVecGa::NormInfinity"
int
TaoVecGa::NormInfinity (TaoScalar * vnorm)
{
  TaoFunctionBegin;
  GA_Norm_infinity (this->pv, vnorm);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Norm1"
int
TaoVecGa::Norm1 (TaoScalar * vnorm)
{
  TaoFunctionBegin;
  GA_Norm1 (this->pv, vnorm);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Norm2"
int
TaoVecGa::Norm2 (TaoScalar * vnorm)
{
  TaoFunctionBegin;
  *vnorm = sqrt (GA_Ddot (this->pv, this->pv));
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Norm2squared"
int
TaoVecGa::Norm2squared (TaoScalar * vnorm)
{
  TaoFunctionBegin;
  *vnorm = GA_Ddot (this->pv, this->pv);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Scale"
int
TaoVecGa::Scale (TaoScalar alpha)
{
  TaoFunctionBegin;
  GA_Scale (this->pv, &alpha);
  TaoFunctionReturn(0);
}


#undef __FUNCT__ 
#define __FUNCT__ "TaoVecGa::Axpy"
int
TaoVecGa::Axpy (TaoScalar alpha, TaoVec * tv)
{
  int info, pv2 = 0;
  TaoScalar a = 1.0, b = alpha;
  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv2); CHKERRQ(info);
  GA_Add (&a, this->pv, &b, pv2, this->pv);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Axpby"
int
TaoVecGa::Axpby (TaoScalar alpha, TaoVec * tv, TaoScalar beta)
{
  int info, pv2 = 0;
  TaoScalar a = alpha, b = beta;
  TaoFunctionBegin;

  info = TaoVecGetGaVec (tv, &pv2); CHKERRQ(info);
  GA_Add (&a, pv2, &b, this->pv, this->pv);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Aypx"
int
TaoVecGa::Aypx (TaoScalar alpha, TaoVec * x)
{
  int info, pv2 = 0;
  TaoScalar a = alpha, b = 1.0;
  TaoFunctionBegin;

  info = TaoVecGetGaVec (x, &pv2); CHKERRQ(info);
  GA_Add (&a, this->pv, &b, pv2, this->pv);
  TaoFunctionReturn(0);
}


#undef __FUNCT__ 
#define __FUNCT__ "TaoVecGa::AddConstant"
int
TaoVecGa::AddConstant (TaoScalar c)
{
  TaoFunctionBegin;
  GA_Add_constant (this->pv, &c);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Dot"
int
TaoVecGa::Dot (TaoVec * tv, TaoScalar * vDotv)
{
  int info, pv2 = 0;
  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv2); CHKERRQ(info);
  *vDotv = GA_Ddot (this->pv, pv2);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Negate"
int
TaoVecGa::Negate ()
{
  TaoScalar m1 = -1.0;
  TaoFunctionBegin;
  GA_Scale (this->pv, &m1);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Reciprocal"
int
TaoVecGa::Reciprocal ()
{
  TaoFunctionBegin;
  GA_Recip (this->pv);
  TaoFunctionReturn(0);
}


#undef __FUNCT__ 
#define __FUNCT__ "TaoVecGa::GetDimension"
int
TaoVecGa::GetDimension (int *n)
{
  int type, ndim, dims;
  TaoFunctionBegin;
  NGA_Inquire (this->pv, &type, &ndim, &dims);
  if (ndim != 1)
    SETERRQ(1,"GA Vector has bad dimension");
  *n = dims;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::PointwiseMultiply"
int
TaoVecGa::PointwiseMultiply (TaoVec * tv, TaoVec * tw)
{
  int info, pv1 = 0, pv2 = 0;
  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv1); CHKERRQ(info);
  info = TaoVecGetGaVec (tw, &pv2); CHKERRQ(info);
  GA_Elem_multiply (pv1, pv2, this->pv);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::PointwiseDivide"
int
TaoVecGa::PointwiseDivide (TaoVec * tv, TaoVec * tw)
{
  int info, pv1 = 0, pv2 = 0;
  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv1); CHKERRQ(info);
  info =  TaoVecGetGaVec (tw, &pv2); CHKERRQ(info);
  GA_Elem_divide (pv1, pv2, this->pv);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Median"
int
TaoVecGa::Median (TaoVec * tv, TaoVec * tw, TaoVec * tx)
{
  int info, pv1 = 0, pv2 = 0, pv3 = 0;
  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv1); CHKERRQ(info);
  info = TaoVecGetGaVec (tw, &pv2); CHKERRQ(info);
  info = TaoVecGetGaVec (tx, &pv3); CHKERRQ(info);

  GA_Median (pv1, pv2, pv3, this->pv);

  TaoFunctionReturn(0);
}

#undef __FUNCT__ 
#define __FUNCT__ "TaoVecGa::PointwiseMinimum"
int
TaoVecGa::PointwiseMinimum (TaoVec * tv, TaoVec * tw)
{
  int info, pv1 = 0, pv2 = 0;
  info = TaoVecGetGaVec (tv, &pv1); CHKERRQ(info);
  info = TaoVecGetGaVec (tw, &pv2); CHKERRQ(info);

  GA_Elem_minimum (pv1, pv2, this->pv);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Fischer"
int
TaoVecGa::Fischer (TaoVec * tx, TaoVec * tf, TaoVec * tl, TaoVec * tu)
{
  TaoFunctionBegin;
  SETERRQ(1," TaoVecGa::Fischer not implemented.");
  TaoFunctionReturn(0);
}

#undef __FUNCT__ 
#define __FUNCT__ "TaoVecGa::PointwiseMaximum"
int
TaoVecGa::PointwiseMaximum (TaoVec * tv, TaoVec * tw)
{
  int info, pv1 = 0, pv2 = 0;
  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv1); CHKERRQ(info);
  info = TaoVecGetGaVec (tw, &pv2); CHKERRQ(info);
  GA_Elem_maximum (pv1, pv2, this->pv);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::Waxpby"
//Sums two scaled vectors and stores the result in this vector. (this=alpha*xx+beta*yy) 
int
TaoVecGa::Waxpby (TaoScalar a, TaoVec * tv, TaoScalar b, TaoVec * tw)
{
  int info, pv1 = 0, pv2 = 0;
  TaoFunctionBegin;
  info = TaoVecGetGaVec (tv, &pv1); CHKERRQ(info);
  info = TaoVecGetGaVec (tw, &pv2); CHKERRQ(info);

  if ((a == 0) && (b == 0))
    GA_Zero (this->pv);
  else if (a == 0)
    {
      GA_Scale (pv2, &b);
      GA_Copy (pv2, this->pv);
    }
  else if (b == 0)
    {
      GA_Scale (pv1, &a);
      GA_Copy (pv1, this->pv);
    }
  else
    GA_Add (&a, pv1, &b, pv2, this->pv);

  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::AbsoluteValue"
int
TaoVecGa::AbsoluteValue ()
{
  TaoFunctionBegin;
  GA_Abs_value (this->pv);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::MinElement"
int
TaoVecGa::MinElement (TaoScalar * val)
{
  int index;
  TaoFunctionBegin;
  NGA_Select_elem (this->pv, (char*)"min", val, &index);
  TaoFunctionReturn(0);
}

/*

int TaoVecGa::GetArray (TaoScalar ** dptr, int *n)
dptr - a pointer to the pointer to the first memory location on the local patch (output)
n    - the number of elements of the GA vector at the location pointed by dptr. (output)
*/
#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::GetArray"
int
TaoVecGa::GetArray (TaoScalar ** dptr, int *n)
{
  int type, ndim, dims;
  int lo, hi, ld;
  TaoScalar *ptr = 0;

  TaoFunctionBegin;
  int me = GA_Nodeid ();

  GA_Sync ();

  NGA_Inquire (this->pv, &type, &ndim, &dims);
  if (type != MT_C_DBL) {
    //the only datatype for taovec is TaoScalar
    SETERRQ(1,"Global Array is wrong data type");
  }

  if (ndim != 1) {
    //we only deal with one dimensional global array as a vector
    SETERRQ(1,"Global Array has wrong dimension");
  }

  NGA_Distribution (this->pv, me, &lo, &hi);

  //Get the size of the patch on this processor
  if (lo >= 0 && hi >= 0) {
    NGA_Access (this->pv, &lo, &hi, &ptr, &ld);
    *n = hi - lo + 1;
  } else 
    *n = 0;

  *dptr = ptr;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::RestoreArray"
int
TaoVecGa::RestoreArray (TaoScalar ** dptr, int *n)
{
  int lo, hi;
  int me = GA_Nodeid ();

  TaoFunctionBegin;
  GA_Sync ();
  NGA_Distribution (this->pv, me, &lo, &hi);

  if (lo != 0 && hi != -1) 
    NGA_Release_update (this->pv, &lo, &hi);

  GA_Sync ();

  *dptr = 0;
  *n = 0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::SetReducedVec"
int
TaoVecGa::SetReducedVec (TaoVec * VV, TaoIndexSet * SS)
{
  TaoFunctionBegin;
  SETERRQ(1,"TaoVecGa::SetReducedVec not yet implemented");
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::ReducedXPY"
int
TaoVecGa::ReducedXPY (TaoVec * VV, TaoIndexSet * SS)
{
  TaoFunctionBegin;
  SETERRQ(1,"TaoVecGa::ReducedXPY not yet implemented.");
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::ReducedCopyFromFull"
int
TaoVecGa::ReducedCopyFromFull (TaoVec * VV, TaoIndexSet * SS)
{
  TaoFunctionBegin;
  SETERRQ(1,"TaoVecGa::ReducedCopyFromFull not yet implemented.");
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::StepMax"
int
TaoVecGa::StepMax (TaoVec * tv, double * step)
{
  int info;
  TaoFunctionBegin;
  GAVec X = this->GetVec ();
  GAVec DX = 0;
  info = TaoVecGetGaVec (tv, &DX); CHKERRQ(info);
  GA_Step_max (X, DX, step);
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::StepMax2"
int
TaoVecGa::StepMax2 (TaoVec * tv, TaoVec * xl, TaoVec * xu, TaoScalar * step)
{

  TaoFunctionBegin;
  SETERRQ(1,"TaoVecGa::CreatIndexSet not yet implemented.");
  TaoFunctionReturn(0);

}



#undef __FUNCT__ 
#define __FUNCT__ "TaoVecGa::View"
int
TaoVecGa::View ()
{
  TaoFunctionBegin;
  GA_Print (this->pv);
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::CreateIndexSet"
int
TaoVecGa::CreateIndexSet (TaoIndexSet ** SSS)
{
  TaoFunctionBegin;
  SETERRQ(1,"TaoVecGa::CreatIndexSet not yet implemented.");
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::StepBoundInfo"
int 
TaoVecGa::StepBoundInfo(TaoVec *TXL, TaoVec *TXU, TaoVec *S, 
			double *bmin1, double *bmin2, double *bmax) {
  int i, n, n2, n3, n4, info;
  double tt,t1=1.0e+20, t2=1.0e+20, t3=0;
  GAScalar *x,*xl,*xu,*dx;

  TaoVecGa *X=this;

  TaoFunctionBegin;
  info = X->GetArray((TaoScalar **)&x, &n); CHKERRQ(info);
  info = TXL->GetArray((TaoScalar **)&xl, &n2); CHKERRQ(info);
  info = TXU->GetArray((TaoScalar **)&xu, &n3); CHKERRQ(info);
  info = S->GetArray((TaoScalar **)&dx, &n4); CHKERRQ(info);

  if ( (n != n2) || (n != n3) || (n != n4) ) {
    SETERRQ(1,"Vectors must all be same size");
  }
  
  for (i=0;i<n;i++){
    if (dx[i]>0){
      tt=(xu[i]-x[i])/dx[i];
      t1=TaoMin(t1,tt);
      if (t1>0){
        t2=TaoMin(t2,tt);
      }
      t3=TaoMax(t3,tt);
    } else if (dx[i]<0){
      tt=(xl[i]-x[i])/dx[i];
      t1=TaoMin(t1,tt);
      if (t1>0){
        t2=TaoMin(t2,tt);
      }
      t3=TaoMax(t3,tt);
    }
  }
  
  info = X->RestoreArray((TaoScalar **)&x, &n); CHKERRQ(info);
  info = TXL->RestoreArray((TaoScalar **)&xl, &n2); CHKERRQ(info);
  info = TXU->RestoreArray((TaoScalar **)&xu, &n3); CHKERRQ(info);
  info = S->RestoreArray((TaoScalar **)&dx, &n4); CHKERRQ(info);
  
  if (bmin1) GA_Dgop(&t1, 1, (char*)"min");
  if (bmin2) GA_Dgop(&t2, 1, (char*)"min");
  if (bmax)  GA_Dgop(&t3, 1, (char*)"max");
  
  *bmin1 = t1; /* boundmin */
  *bmin2 = t2; /* wolfemin */
  *bmax  = t3; /* boundmax */
  
  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "TaoVecGa::BoundGradientProjection"
int 
TaoVecGa::BoundGradientProjection(TaoVec *GG, TaoVec *XXL, TaoVec *XX, 
				  TaoVec *XXU) {
  
  TaoVecGa *GP  = this;
  
  int i, n, n2, n3, n4, n5, info; 
  GAScalar *xptr,*xlptr,*xuptr,*gptr,*gpptr;
  GAScalar xval,gpval;
  
  TaoFunctionBegin;
  /* Project variables at the lower and upper bound */   

  info = XX->GetArray((TaoScalar **)&xptr, &n); CHKERRQ(info);
  info = XXL->GetArray((TaoScalar **)&xlptr, &n2); CHKERRQ(info);
  info = XXU->GetArray((TaoScalar **)&xuptr, &n3); CHKERRQ(info);
  info = GG->GetArray((TaoScalar **)&gptr, &n4); CHKERRQ(info);

  if ( (n != n2) || (n != n3) || (n != n4) ) {
    SETERRQ(1,"Vectors must all be same size");
  }
  

  if (GG!=GP){
    info = GP->GetArray((TaoScalar **)&gpptr,&n5);  CHKERRQ(info);
    if (n != n5) {SETERRQ(1, "Vectors must all be same size");}
  } else { gpptr=gptr; }

  for (i=0; i<n; ++i){
    gpval = gptr[i]; xval = xptr[i];
    
    if (gpval>0 && xval<=xlptr[i]){
      gpval = 0;
    } else if (gpval<0 && xval>=xuptr[i]){
      gpval = 0;
    }
    gpptr[i] = gpval;
  }
  
  info = XX->RestoreArray((TaoScalar **)&xptr, &n); CHKERRQ(info);
  info = XXL->RestoreArray((TaoScalar **)&xlptr, &n2); CHKERRQ(info);
  info = XXU->RestoreArray((TaoScalar **)&xuptr, &n3); CHKERRQ(info);
  info = GG->RestoreArray((TaoScalar **)&gptr, &n4); CHKERRQ(info);
  if (GG!=GP) {
    info = GP->RestoreArray((TaoScalar **)&gpptr, &n5); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}


