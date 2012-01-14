//File: taomat_ga.c

/**************************************************************

Author: Limin Zhang, Ph.D.
        Mathematics Department
        Columbia Basin College
        Pasco, WA 99301
        Limin.Zhang@cbc2.org

Mentor: Jarek Naplocha, Ph.D.
        Environmental Molecular Science Laboratory
        Pacific Northwest National Laboratory
        Richland, WA 99352

Date: 7/11/2002

Purpose:
      to design and implement some interfaces between TAO and
      global arrays.

Revise History:

8/8/02
	To clean Petsc function calls and marcos.
**************************************************************/



#include "taomat_ga.h"
#include "taovec_ga.h"
#include <string>
/*
int MatD_Fischer(Mat, Vec, Vec, Vec, Vec, Vec, Vec, Vec, Vec);
int MatCreateADA(Mat,Vec,Vec,Mat*);
int MatSMFResetRowColumn(Mat,IS,IS);
*/

/*@C
   TaoWrapGaMat - Creates a new TaoMat object using a GA matrix.

   Input Parameter:
.  M -  a GA matrix

   Output Parameter:
.  MM - new TaoMat

   Level: advanced

.seealso TaoMatGetGaMat(), TaoMatDestroy()
@*/
int TaoWrapGaMat (GAMat M, TaoMatGa ** MM)
{
  TaoFunctionBegin;
  *MM = new TaoMatGa (M);
  TaoFunctionReturn (0);
}

/*@C
   TaoMatGetGaMat - If the TaoMat is of the TaoMatGa type, this routine
   gets the underlying GA matrix.

   Input Parameter:
.  MM - the TaoMat 

   Output Parameter:
.  M -  the GA mat

   Note:
   The function TaoMatGa::GetMat() will also return the Mat M.

   Level: advanced
@*/
int TaoMatGetGaMat (TaoMat * MM, GAMat * M)
{
  TaoFunctionBegin;
  if (MM)
    {
      *M = ((TaoMatGa *) MM)->GetMat ();
    }
  else
    {
      *M = 0;
    }
  TaoFunctionReturn (0);
}

TaoMatGa::TaoMatGa (GAMat MM):TaoMat ()
{
  this->pm = MM;
  //  this->MatObject = (void *) (&MM);
  this->pm_pre = MM;
  return;
}

int TaoMatGa::Compatible (TaoVec * xx, TaoVec * yy, TaoTruth *flag)
{
  int
    nx,
    ny,
    m,
    n;
  GAVec
    x,
    y;
  if (xx == 0 || yy == 0)
    { *flag=TAO_FALSE; return 0; }

  x = ((TaoVecGa *) xx)->GetVec ();
  y = ((TaoVecGa *) yy)->GetVec ();
  if (x == 0 || y == 0 || this->pm == 0)
    { *flag=TAO_FALSE; return 0; }

  int
    xtype,
    ytype,
    type,
    ndim,
    dims[2];

  NGA_Inquire (x, &xtype, &ndim, &nx);
  if (ndim != 1)
    GA_Error ((char*)"TaoMatGa::Compatible: Wrong dimension", ndim);

  NGA_Inquire (y, &ytype, &ndim, &ny);
  if (ndim != 1)
    GA_Error ((char*)"TaoMatGa::Compatible:Wrong dimension", ndim);

  NGA_Inquire (this->pm, &type, &ndim, dims);
  if (ndim != 2)
    GA_Error ((char*)"TaoMatGa::Compatible:Wrong dimension", ndim);
  n = dims[0];
  m = dims[1];

  if (xtype != type || ytype != type)
    GA_Error ((char*)"TaoMatGa::Compatible:Wrong data type", type);

  if (n != nx || m != ny)
    { *flag=TAO_FALSE; return 0; }

  *flag=TAO_TRUE; 
  return 0;

}

int
TaoMatGa::Clone (TaoMat ** ntm)
{
  TaoMatGa *nptm;
  GAMat npm;
  TaoFunctionBegin;
  npm = GA_Duplicate (this->pm, (char*)"Cloned");	
  if (!npm)
    GA_Error ((char*)"TaoMatGa::Clone:duplicate failed: ", npm);
  TaoWrapGaMat (npm, &nptm);
  *ntm = nptm;
  TaoFunctionReturn (0);
}

int
TaoMatGa::CopyFrom (TaoMat * tm)
{
  TaoFunctionBegin;
#if 0 //try a different way to get GAMat
  TaoMatGa *temp2 = dynamic_cast < TaoMatGa * >(tm);
  GAMat pm2 = temp2->GetMat ();
#endif
  GAMat pm2 = ((TaoMatGa *)tm)->GetMat ();
  GA_Copy (pm2, this->pm);
  TaoFunctionReturn (0);
}

int
TaoMatGa::GetDimensions (int *m, int *n)
{
  TaoFunctionBegin;
  int type, ndim, dims[2];
  NGA_Inquire (this->pm, &type, &ndim, dims);
  if (ndim != 2)
    GA_Error ((char*)"TaoMatGa::GetDimension: wrong dimension", ndim);
  else
    {
      *n = dims[0];
      *m = dims[1];
    }
  TaoFunctionReturn (0);
}

int
TaoMatGa::Multiply (TaoVec * tv, TaoVec * tw)
{
  double alpha = 1.0;
  double beta = 0.0;
  TaoTruth flag;
  TaoFunctionBegin;
  GAVec vv = ((TaoVecGa *) tv)->GetVec ();
  GAVec ww = ((TaoVecGa *) tw)->GetVec ();
  int m, n = 1, k;		//n = 1 due to vecotor ww that is k by 1 array 
  this->GetDimensions (&m, &k);
  this->Compatible (tv, tw,&flag);
  if (flag)
    GA_Dgemm ('N', 'N', m, n, k, alpha, this->pm, vv, beta, ww);
  else
    TaoFunctionReturn (1);
  TaoFunctionReturn (0);
}

int
TaoMatGa::MultiplyAdd (TaoVec * tv, TaoVec * tw, TaoVec * ty)
{
  TaoTruth flag1,flag2;
  TaoFunctionBegin;
  GAVec xx = ((TaoVecGa *) tv)->GetVec ();
  GAVec ww = ((TaoVecGa *) tw)->GetVec ();
  GAVec yy = ((TaoVecGa *) ty)->GetVec ();
  GAVec temp_ww = GA_Duplicate (ww, (char*)"temp_WW");
  if (!temp_ww)
    GA_Error ((char*)"TaoMatGa::MultiplyAdd:Duplicate:", temp_ww);
  GA_Copy (ww, temp_ww);
  double alpha = 1.0;
  double beta = 1.0;
  int m, n = 1, k;		//n = 1 due to vecotor ww that is k by 1 array
  this->GetDimensions (&m, &k);
  this->Compatible (tv, ty, &flag1);
  this->Compatible (tv, tw, &flag2);
  if (flag1 && flag2)
    {
      GA_Dgemm ('N', 'N', m, n, k, alpha, this->pm, xx, beta, temp_ww);
      //copy from temp_ww to yy. Doing this way, the data in ww is not dirty!!!!!
      GA_Copy (temp_ww, yy);
      //free the temporary memory
      GA_Destroy(temp_ww);
    }
  else
    TaoFunctionReturn (1);
  TaoFunctionReturn (0);
}

int
TaoMatGa::MultiplyTranspose (TaoVec * tv, TaoVec * tw)
{
  double alpha = 1.0;
  double beta = 0.0;
  int m, n = 1, k;		/*n = 1 due to vecotor ww that is k */
  TaoFunctionBegin;
  GAVec vv = ((TaoVecGa *) tv)->GetVec ();
  GAVec ww = ((TaoVecGa *) tw)->GetVec ();
  this->GetDimensions (&m, &k);
  GA_Dgemm ('T', 'N', k, n, m, alpha, this->pm, vv, beta, ww);
  TaoFunctionReturn (0);
}

int
TaoMatGa::MultiplyTransposeAdd (TaoVec * tv, TaoVec * tw, TaoVec * ty)
{
  GAVec xx = ((TaoVecGa *) tv)->GetVec ();
  GAVec ww = ((TaoVecGa *) tw)->GetVec ();
  GAVec yy = ((TaoVecGa *) ty)->GetVec ();
  TaoFunctionBegin;
  TaoFunctionBegin;
  GAVec temp_ww = GA_Duplicate (ww, (char*)"temp_WW");
  if (!temp_ww)
    GA_Error ((char*)"TaoMatGa::MultiplyAdd:Dupli cate:", temp_ww);
  GA_Copy (ww, temp_ww);
  double alpha = 1.0;
  double beta = 1.0;
  int m, n = 1, k;
  //n = 1 due to vecotor ww that is k by 1 array
  this->GetDimensions (&m, &k);
  GA_Dgemm ('T', 'N', k, n, m, alpha, this->pm, xx, beta, temp_ww);
  //copy from temp_ww to yy. Doing this way, the data in ww is not dirty!!!!!
  GA_Copy (temp_ww, yy);
  //free the memory held by the temporary GA array temp_ww
  GA_Destroy(temp_ww);
  TaoFunctionReturn (0);
}

int
TaoMatGa::SetDiagonal (TaoVec * tv)
{
  GAVec vv = ((TaoVecGa *) tv)->GetVec ();
  TaoFunctionBegin;
  GA_Set_diagonal (this->pm, vv);
  TaoFunctionReturn (0);
}

int
TaoMatGa::AddDiagonal (TaoVec * tv)
{
  GAVec vv = ((TaoVecGa *) tv)->GetVec ();
  TaoFunctionBegin;
  GA_Add_diagonal (this->pm, vv);
  TaoFunctionReturn (0);
}

int
TaoMatGa::GetDiagonal (TaoVec * tv)
{
  GAVec vv = ((TaoVecGa *) tv)->GetVec ();
  TaoFunctionBegin;
  GA_Get_diag (this->pm, vv);
  TaoFunctionReturn (0);
}

int
TaoMatGa::ShiftDiagonal (TaoScalar c)
{
  TaoFunctionBegin;
  GA_Shift_diagonal (this->pm, &c);
  TaoFunctionReturn (0);
}

int
TaoMatGa::View ()
{
  TaoFunctionBegin;
  GA_Print (this->pm);
  TaoFunctionReturn (0);
}


int
TaoMatGa::Presolve ()
{
  TaoFunctionBegin;
  TaoFunctionReturn (0);
}

int
TaoMatGa::Solve (TaoVec * tv, TaoVec * tw, TaoTruth * tt)
{
  TaoFunctionBegin;
  GA_Error ((char*)"not yet implemented.", 0);
  TaoFunctionReturn (1);
}

int
TaoMatGa::RowScale (TaoVec * tv)
{
  GAVec vv = ((TaoVecGa *) tv)->GetVec ();
  TaoFunctionBegin;
  GA_Scale_rows (this->pm, vv);
  TaoFunctionReturn (0);
}

int
TaoMatGa::ColScale (TaoVec * tv)
{
  GAVec vv = ((TaoVecGa *) tv)->GetVec ();
  TaoFunctionBegin;
  GA_Scale_cols (this->pm, vv);
  TaoFunctionReturn (0);
}

int
TaoMatGa::CreateReducedMatrix (TaoIndexSet * S1,
			       TaoIndexSet * S2, TaoMat ** MM)
{
  int info;
  TaoMatGa *M;
  GAMat B = 0;
  TaoFunctionBegin;
  info = TaoWrapGaMat (B, &M);
  CHKERRQ (info);
  *MM = M;
  TaoFunctionReturn (0);
}

int
TaoMatGa::SetReducedMatrix (TaoMat * M, TaoIndexSet * S1, TaoIndexSet * S2)
{
  TaoFunctionBegin;
  GA_Error ((char*)"Not yet implemented.", 0);
  TaoFunctionReturn (1);
}



int
TaoMatGa::D_Fischer (TaoVec * tx, TaoVec * tf, TaoVec * tl,
		     TaoVec * tu, TaoVec * tda,
		     TaoVec * tdb, TaoVec * tt1, TaoVec * tt2)
{
  TaoFunctionBegin;
  GA_Error ((char*)"Not yet implemented.", 0);
  TaoFunctionReturn (1);
}

int
TaoMatGa::Norm1 (double *nm)
{
  TaoFunctionBegin;
  GA_Norm1 (this->pm, nm);
  TaoFunctionReturn (0);
}


int
TaoMatGa::Fill (double c)
{

  TaoFunctionBegin;
  GA_Fill (this->pm, &c);
  TaoFunctionReturn (0);
}


int
TaoMatGa::Zero ()
{

  TaoFunctionBegin;
  GA_Zero (this->pm);
  TaoFunctionReturn (0);

}



