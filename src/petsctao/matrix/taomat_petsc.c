#include "tao_general.h"

#ifdef TAO_USE_PETSC

#include "taomat_petsc.h"
#include "../vector/taovec_petsc.h"
#include "../indexset/taois_petsc.h"

int MatD_Fischer(Mat, Vec, Vec, Vec, Vec, Vec, Vec, Vec, Vec);
int MatD_SFischer(Mat, Vec, Vec, Vec, Vec, double, Vec, Vec, Vec, Vec, Vec);
int MatSMFResetRowColumn(Mat,IS,IS);

#undef __FUNCT__
#define __FUNCT__ "TaoWrapPetscMat"
/*@C
   TaoWrapPetscMat - Creates a new TaoMat object using a PETSc matrix.

   Input Parameter:
+  M -  a PETSc matrix
-  MM -  the address of a pointer to a TaoMatPetsc

   Output Parameter:
.  MM - the address of a pointer to new TaoMat

   Note:  
   A TaoMatPetsc is an object with the methods of an abstract
   TaoMat object.  A TaoMatPetsc contains an implementation of the TaoMat
   methods.  Routines using these vectors should declare a pointer to 
   a TaoMat, assign this pointer to the address of a TaoMat object, 
   use the pointer to invoke methods on the object, and use this pointer
   as an argument when calling other routines.  This usage is different
   from the usage of a PETSc Mat.  In PETSc, applications will typically
   declare a Mat, and pass it as an argument into routines.  That is,
   applications will typically declare a pointer to a TaoMat and use the
   pointer, or declare a Mat and use it directly.

   Note:
   The user is repsonsible for destroying the Mat M, in addition to
   to TaoMatPetsc vector MM.  The Mat can be destroyed immediately
   after this routine.

  Level: developer

.seealso TaoMatGetPetscMat(), TaoMatDestroy()
@*/
int TaoWrapPetscMat( Mat M, TaoMatPetsc* *MM){
  PetscFunctionBegin;
  if (MM){ *MM = new  TaoMatPetsc(M);}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::TaoMatPetsc"
/* @C
   TaoMatPetsc - Creates a new TaoMat object using a PETSc matrix.

   Input Parameter:
+  MM -  a PETSc matrix

   Level: advanced

   Note:
   The method TaoMatPetsc::SetMatrix(Mat M, Mat MPre, MatStructure) replaces the matrix MM with M.
   The method TaoMatPetsc::GetMatrix() returns a pointer to the matrix MM, its preconditioner, and KSP flag

@ */
TaoMatPetsc::TaoMatPetsc(Mat MM) : TaoMat()
{
  int info;
  pm=0; pm_pre=0;
  preflag=DIFFERENT_NONZERO_PATTERN;
  MPI_Comm comm;
  if (MM) { 
    PetscObjectGetComm((PetscObject)MM,&comm);
    pmatviewer=PETSC_VIEWER_STDOUT_(comm);
    SetMatrix(MM,MM,preflag);  
    info = PetscInfo(MM,"Wrap a PETSc Mat to create a TaoMat .\n");
  }
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::~TaoMatPetsc"
TaoMatPetsc::~TaoMatPetsc()
{
  int info;
  if (pm) {
    info = PetscInfo(pm,"TAO: Destroy a PETSc Mat within a TaoMat .\n");
    PetscObjectDereference((PetscObject)pm);
    PetscObjectDereference((PetscObject)pm_pre);
  }
  pm=0; pm_pre=0;
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::GetMatrix"
/* @C
   GetMatrix - Gets the underlying PETSc matrix information

   Output Parameters:
+  M -  a PETSc matrix
.  MPre -  a PETSc preconditioning matrix (use PETSC_NULL for no change)
-  flag -  flag associated with KSPSetOperators (PETSC_NULL for no change).

   Level: advanced

.seealso TaoMatPetsc::SetMatrix()
 
@ */
int TaoMatPetsc::GetMatrix(Mat *M, Mat *MPre, MatStructure *flag){
  PetscFunctionBegin;
  if (M) *M = pm;
  if (flag) *flag = preflag;
  if (MPre) *MPre = pm_pre;
  PetscFunctionReturn(0);  
}


#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::SetMatrix"
/* @C
   SetMatrix - Changes the underlying PETSc matrix structrure

   Input Parameters:
+  M -  a PETSc matrix
.  MPre -  a PETSc preconditioning matrix (use PETSC_NULL for no change)
-  flag -  flag associated with KSPSetOperators (PETSC_NULL for no change).

   Level: advanced

   Note:
   The method TaoMatPetsc::GetMatrix() returns a pointer to the matrix M

@ */
int TaoMatPetsc::SetMatrix(Mat M, Mat MPre, MatStructure flag){
  int info;
  PetscInt r1,r2,c1,c2,R1,R2,C1,C2;
  PetscFunctionBegin;
  if (M || MPre){
    PetscValidHeaderSpecific(M,MAT_COOKIE,1);
    PetscValidHeaderSpecific(MPre,MAT_COOKIE,2);
    if (M!=MPre){
      info = MatGetSize(M,&R1,&C1);CHKERRQ(info);
      info = MatGetSize(MPre,&R2, &C2);CHKERRQ(info);
      info = MatGetLocalSize(M,&r1,&c1);CHKERRQ(info);
      info = MatGetLocalSize(MPre,&r2,&c2);CHKERRQ(info);
      if (R1!=R2 || C1!=C2 || r1!=r2 || c1!=c2) {
	SETERRQ(1,"TAO Error: Preconditioning Matrix does not match the Solution matrix");
      }
    }
    
    info = PetscObjectReference((PetscObject)M);CHKERRQ(info);
    info = PetscObjectReference((PetscObject)MPre);CHKERRQ(info); 
  }

  if (pm&&pm_pre){
    info = PetscObjectDereference((PetscObject)pm);CHKERRQ(info);
    info = PetscObjectDereference((PetscObject)pm_pre);CHKERRQ(info); 
  }

  pm=M;
  pm_pre=MPre; 
  preflag=flag;
  PetscFunctionReturn(0);  
}



#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::Compatible"
int TaoMatPetsc::Compatible(TaoVec* xx, TaoVec* yy, TaoTruth *flag){
  int info;
  PetscInt pmCol,pmRow,pmrow,pmcol,xN,yN,yn;
  Vec x,y;

  PetscFunctionBegin;
  *flag=TAO_TRUE;
  if (xx==0 || yy==0){
    *flag=TAO_FALSE;
    PetscFunctionReturn(0);
  }
  x=((TaoVecPetsc*)xx)->GetVec();
  y=((TaoVecPetsc*)yy)->GetVec();
  info = MatGetSize(pm,&pmRow,&pmCol); CHKERRQ(info);
  info = MatGetLocalSize(pm,&pmrow,&pmcol); CHKERRQ(info);
  info = VecGetSize(x,&xN); CHKERRQ(info);
  info = VecGetSize(y,&yN); CHKERRQ(info);
  info = VecGetLocalSize(y,&yn); CHKERRQ(info);
  if ((pmCol != xN) || (pmRow != yN) || (pmrow != yn)) {
    *flag=TAO_FALSE;
    PetscFunctionReturn(0);
  }
  
  PetscValidHeaderSpecific(x,VEC_COOKIE,1); 
  PetscValidHeaderSpecific(y,VEC_COOKIE,2); 
  PetscValidHeaderSpecific(pm,MAT_COOKIE,0); 
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::Clone"
int TaoMatPetsc::Clone(TaoMat* *ntm){
  int info;
  TaoMatPetsc *nptm;
  Mat M,MPre,npm,nmpre;
  MatStructure flag;

  PetscFunctionBegin;
  info = GetMatrix(&M,&MPre,&flag); CHKERRQ(info);
  info = MatDuplicate(M,MAT_COPY_VALUES,&npm); CHKERRQ(info);
  info = MatAssemblyBegin(npm,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(npm,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  if (M!=MPre){
    info = MatDuplicate(MPre,MAT_COPY_VALUES,&nmpre); CHKERRQ(info);
    info = MatAssemblyBegin(nmpre,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
    info = MatAssemblyEnd(nmpre,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  }
  info = TaoWrapPetscMat(npm, &nptm);   CHKERRQ(info);
  *ntm = nptm;
  info = nptm->SetMatrix(npm,nmpre,flag); CHKERRQ(info);
  info = MatDestroy(npm); CHKERRQ(info);
  if (M!=MPre){
    info = MatDestroy(nmpre); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::CopyFrom"
int TaoMatPetsc::CopyFrom(TaoMat* tm){
  int info;
  TaoMatPetsc* tmm=(TaoMatPetsc*)tm;
  Mat M,MPre;
  MatStructure flag;

  PetscFunctionBegin;
  info = tmm->GetMatrix(&M,&MPre,&flag);CHKERRQ(info);
  info = MatCopy(M,pm,SAME_NONZERO_PATTERN); 
  CHKERRQ(info);
  if (pm!=pm_pre){
    info = MatCopy(MPre,pm_pre,SAME_NONZERO_PATTERN); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::GetDimensions"
int TaoMatPetsc::GetDimensions( TaoInt *m, TaoInt *n ){
  int info;
  PetscFunctionBegin;
  info=MatGetSize(pm,m,n); CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::Multiply"
int TaoMatPetsc::Multiply(TaoVec *tv, TaoVec *tw)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *pw = dynamic_cast <TaoVecPetsc *> (tw);

  int info;
  PetscFunctionBegin;
  info=MatMult(pm,pv->GetVec(),pw->GetVec()); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::MultiplyAdd"
int TaoMatPetsc::MultiplyAdd(TaoVec* tv,TaoVec* tw,TaoVec* ty)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *pw = dynamic_cast <TaoVecPetsc *> (tw);
  TaoVecPetsc *py = dynamic_cast <TaoVecPetsc *> (ty);
  int info;
  PetscFunctionBegin;
  info=MatMultAdd(pm,pv->GetVec(),pw->GetVec(),py->GetVec()); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::MultiplyTranspose"
int TaoMatPetsc::MultiplyTranspose(TaoVec* tv,TaoVec* tw)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *pw = dynamic_cast <TaoVecPetsc *> (tw);
  int info;
  PetscFunctionBegin;
  info=MatMultTranspose(pm,pv->GetVec(),pw->GetVec()); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::MultiplyTransposeAdd"
int TaoMatPetsc::MultiplyTransposeAdd(TaoVec* tv,TaoVec* tw,TaoVec* ty)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  TaoVecPetsc *pw = dynamic_cast <TaoVecPetsc *> (tw);
  TaoVecPetsc *py = dynamic_cast <TaoVecPetsc *> (ty);
  int info;
  PetscFunctionBegin;
  info=MatMultTransposeAdd(pm,pv->GetVec(),pw->GetVec(),py->GetVec()); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::SetDiagonal"
int TaoMatPetsc::SetDiagonal(TaoVec* tv)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  int info;
  PetscFunctionBegin;
  info = MatDiagonalSet(pm,pv->GetVec(),INSERT_VALUES); CHKERRQ(info);
  if (pm != pm_pre){
    info = MatDiagonalSet(pm_pre,pv->GetVec(),INSERT_VALUES); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::AddDiagonal"
int TaoMatPetsc::AddDiagonal(TaoVec* tv)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  int info;
  PetscFunctionBegin;
  info = MatDiagonalSet(pm,pv->GetVec(),ADD_VALUES); CHKERRQ(info);
  if (pm != pm_pre){
    info = MatDiagonalSet(pm_pre,pv->GetVec(),ADD_VALUES); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::GetDiagonal"
int TaoMatPetsc::GetDiagonal(TaoVec* tv)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  int info;
  PetscFunctionBegin;
  info = MatGetDiagonal(pm,pv->GetVec()); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::ShiftDiagonal"
int TaoMatPetsc::ShiftDiagonal(double c)
{
  int info;
  PetscScalar cc=c;
  PetscFunctionBegin;
  info = MatShift(pm,cc); CHKERRQ(info);
  if (pm != pm_pre) {
    info = MatShift(pm_pre,cc); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::SetPetscViewer"
/*@C
   SetPetscViewer

   Input Parameter:
.  viewer - a viewer object to be used with View() and MatView()

   Level: advanced

@*/
int TaoMatPetsc::SetPetscViewer(PetscViewer viewer){
  PetscFunctionBegin;  
  pmatviewer= viewer;
  PetscFunctionReturn(0);  
}


#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::View"
int TaoMatPetsc::View()
{
  int info;
  PetscTruth flg=PETSC_FALSE;
  PetscFunctionBegin;
  info = PetscOptionsHasName(PETSC_NULL,"-tao_mat_draw",&flg); CHKERRQ(info);
  if (flg){
    info = MatView(pm,PETSC_VIEWER_DRAW_WORLD); CHKERRQ(info);
  } else {
    info = MatView(pm,pmatviewer); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::RowScale"
int TaoMatPetsc::RowScale(TaoVec* tv)
{
  TaoVecPetsc *pv = dynamic_cast <TaoVecPetsc *> (tv);
  int info;

  PetscFunctionBegin;
  info = MatDiagonalScale(pm, pv->GetVec(), PETSC_NULL); CHKERRQ(info);
  if (pm!=pm_pre){
    info = MatDiagonalScale(pm_pre,pv->GetVec(),PETSC_NULL); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::ColScale"
int TaoMatPetsc::ColScale(TaoVec* tv){
  Vec vv = ((TaoVecPetsc *)tv)->GetVec();
  int info;

  PetscFunctionBegin;
  info = MatDiagonalScale(pm, PETSC_NULL, vv); CHKERRQ(info);
  if (pm!=pm_pre){
    info = MatDiagonalScale(pm_pre,PETSC_NULL,vv); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::CreateReducedMatrix"
int TaoMatPetsc::CreateReducedMatrix(TaoIndexSet* S1,TaoIndexSet* S2,TaoMat* *MM){
  TaoIndexSetPetsc *prows = dynamic_cast <TaoIndexSetPetsc *> (S1);
  TaoIndexSetPetsc *pcols = dynamic_cast <TaoIndexSetPetsc *> (S2);
  int info;
  TaoMatPetsc *M;
  Mat B=0,BPre,A,APre;
  TaoPetscISType type;
  MatStructure flag;
  PetscTruth assembled,flg;

  PetscFunctionBegin;
  info = prows->GetReducedType(&type); CHKERRQ(info);
  if (type==TaoMaskFullSpace){
    info = GetMatrix(&B,&BPre,&flag); CHKERRQ(info);
    PetscOptionsHasName(0,"-different_submatrix",&flg);
    if (flg==PETSC_TRUE){
      info = TaoWrapPetscMat(B, &M); CHKERRQ(info);
      info = MatAssembled(B,&assembled); CHKERRQ(info);
      if (assembled==PETSC_TRUE){
	info = MatDuplicate(B,MAT_DO_NOT_COPY_VALUES,&A); CHKERRQ(info);
	info = M->SetMatrix(A,A,flag); CHKERRQ(info);
	info = MatDestroy(A); CHKERRQ(info);
      }
      info = MatAssembled(BPre,&assembled); CHKERRQ(info);
      if (B != BPre && assembled==PETSC_TRUE){
	info = MatDuplicate(BPre,MAT_DO_NOT_COPY_VALUES,&APre); CHKERRQ(info);
	info = M->SetMatrix(A,APre,flag); CHKERRQ(info);
	info = MatDestroy(APre); CHKERRQ(info);
      }
    } else {
      info = GetMatrix(&B,&BPre,&flag); CHKERRQ(info);
      info = TaoWrapPetscMat(B, &M); CHKERRQ(info);
      info = M->SetMatrix(B,BPre,flag); CHKERRQ(info);
    }    
  } else if (type==TaoMatrixFree){
    info = GetMatrix(&B,&BPre,&flag); CHKERRQ(info);
    info = MatCreateSubMatrixFree(B,prows->GetIS(), pcols->GetIS(), &A); CHKERRQ(info);
    info = TaoWrapPetscMat(A, &M); CHKERRQ(info);
    info = M->SetMatrix(A,BPre,flag); CHKERRQ(info);
    info = PetscObjectDereference((PetscObject)A);CHKERRQ(info);
  } else {
    info = GetMatrix(&B,&BPre,&flag); CHKERRQ(info);
    info = TaoWrapPetscMat(B, &M); CHKERRQ(info);
    info = M->SetMatrix(B,BPre,flag); CHKERRQ(info);
  }
  *MM=M;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::SetReducedMatrix"
int TaoMatPetsc::SetReducedMatrix(TaoMat* M,TaoIndexSet* S1,TaoIndexSet* S2)
{
  TaoIndexSetPetsc *prows = dynamic_cast <TaoIndexSetPetsc *> (S1);
  TaoIndexSetPetsc *pcols = dynamic_cast <TaoIndexSetPetsc *> (S2);

  int info;
  TaoPetscISType type;
  TaoMatPetsc *MM=((TaoMatPetsc *)M);
  Mat A,Apre,B,BPre;
  MatStructure flag;

  PetscFunctionBegin;
  info=GetMatrix(&B,&BPre,0); CHKERRQ(info);
  info=MM->GetMatrix(&A,&Apre,0); CHKERRQ(info);
  info = prows->GetReducedType(&type); CHKERRQ(info);
  if (type==TaoMaskFullSpace){
    Vec DDiag,VRow,VCol;
    info = GetMatrix(&B,&BPre,&flag); CHKERRQ(info);
    info = MM->GetMatrix(&A,&Apre,&flag); CHKERRQ(info);
    if (!A){
      info = MatDuplicate(B,MAT_COPY_VALUES,&A); CHKERRQ(info);    
      if (B!=BPre){info = MatDuplicate(BPre,MAT_COPY_VALUES,&Apre); CHKERRQ(info);}
      info = MM->SetMatrix(A,Apre,flag); CHKERRQ(info);
      info=MatDestroy(A); CHKERRQ(info);
      if (B!=BPre){info=MatDestroy(Apre); CHKERRQ(info);}
    }
    if (A!=B){ info=MatCopy(A,B,SAME_NONZERO_PATTERN); CHKERRQ(info); }
    if (B!=BPre && Apre!=BPre){ info=MatCopy(Apre,BPre,SAME_NONZERO_PATTERN); CHKERRQ(info); }
    info=prows->GetMask(&VRow); CHKERRQ(info);
    info=pcols->GetMask(&VCol); CHKERRQ(info);
    info=VecDuplicate(VRow,&DDiag); CHKERRQ(info);
    info=MatGetDiagonal(A,DDiag); CHKERRQ(info);
    info=MatDiagonalScale(B,VRow,VCol); CHKERRQ(info);
    info=MatDiagonalSet(B,DDiag,INSERT_VALUES); CHKERRQ(info);
    if (B!=BPre){
      info=MatGetDiagonal(Apre,DDiag); CHKERRQ(info);
      info=MatDiagonalScale(BPre,VRow,VCol); CHKERRQ(info);
      info=MatDiagonalSet(BPre,DDiag,INSERT_VALUES); CHKERRQ(info);
    }
    info=VecDestroy(DDiag); CHKERRQ(info);
  } else if (type==TaoMatrixFree){
    info = MM->GetMatrix(&A,&Apre,&flag); CHKERRQ(info);
    info = MatCreateSubMatrixFree(A,prows->GetIS(),pcols->GetIS(),&B); CHKERRQ(info);
    info = SetMatrix(B,Apre,flag); CHKERRQ(info);
    info = PetscObjectDereference((PetscObject)B);CHKERRQ(info);
  } else {
    IS LocalRows, LocalCols, AllCols;
    PetscInt nn;
    //    info = GetMatrix(&B,&BPre,&flag); CHKERRQ(info);
    info = MM->GetMatrix(&A,&Apre,&flag); CHKERRQ(info);
    info = prows->RedistributeIS(&LocalRows); CHKERRQ(info);
    info = pcols->RedistributeIS(&LocalCols); CHKERRQ(info);
    info = pcols->GetWholeIS(&AllCols); CHKERRQ(info);
    info = ISGetLocalSize(LocalCols,&nn); CHKERRQ(info);
    info = MatGetSubMatrix(A,LocalRows,AllCols,MAT_INITIAL_MATRIX,&B);
    CHKERRQ(info);
    BPre=B;

    if (A != Apre){      
      info = MatGetSubMatrix(Apre,LocalRows,AllCols,MAT_INITIAL_MATRIX,&BPre);
      CHKERRQ(info);
    } else {
      info = PetscObjectReference((PetscObject)BPre);CHKERRQ(info);
    }
    
    info = SetMatrix(B,BPre,flag); CHKERRQ(info);
    info = PetscObjectDereference((PetscObject)B);CHKERRQ(info);
    info = PetscObjectDereference((PetscObject)BPre);CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::D_Fischer"
int TaoMatPetsc::D_Fischer(TaoVec* tx, TaoVec* tf, TaoVec* tl, TaoVec* tu,
                           TaoVec* tt1, TaoVec* tt2, TaoVec* tda, TaoVec* tdb) 
{
  Vec x = ((TaoVecPetsc *)tx)->GetVec();
  Vec f = ((TaoVecPetsc *)tf)->GetVec();
  Vec l = ((TaoVecPetsc *)tl)->GetVec();
  Vec u = ((TaoVecPetsc *)tu)->GetVec();
  Vec da = ((TaoVecPetsc *)tda)->GetVec();
  Vec db = ((TaoVecPetsc *)tdb)->GetVec();
  Vec t1 = ((TaoVecPetsc *)tt1)->GetVec();
  Vec t2 = ((TaoVecPetsc *)tt2)->GetVec();
  int info;

  PetscFunctionBegin;
  info = MatD_Fischer(pm, x, f, l, u, t1, t2, da, db); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::D_SFischer"
int TaoMatPetsc::D_SFischer(TaoVec* tx, TaoVec* tf, TaoVec* tl, TaoVec* tu,
                            double mu,
			    TaoVec* tt1, TaoVec* tt2,
                            TaoVec* tda, TaoVec* tdb, TaoVec* tdm)
{
  Vec x = ((TaoVecPetsc *)tx)->GetVec();
  Vec f = ((TaoVecPetsc *)tf)->GetVec();
  Vec l = ((TaoVecPetsc *)tl)->GetVec();
  Vec u = ((TaoVecPetsc *)tu)->GetVec();
  Vec t1 = ((TaoVecPetsc *)tt1)->GetVec();
  Vec t2 = ((TaoVecPetsc *)tt2)->GetVec();
  Vec da = ((TaoVecPetsc *)tda)->GetVec();
  Vec db = ((TaoVecPetsc *)tdb)->GetVec();
  Vec dm = ((TaoVecPetsc *)tdm)->GetVec();
  int info;

  PetscFunctionBegin;
  if ((mu >= -TAO_EPSILON) && (mu <= TAO_EPSILON)) {
    tdm->SetToZero();
    D_Fischer(tx, tf, tl, tu, tt1, tt2, tda, tdb);
  }
  else {
    info = MatD_SFischer(pm, x, f, l, u, mu, t1, t2, da, db, dm); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMatPetsc::Norm1"
int TaoMatPetsc::Norm1(double *nm){
  int info;
  PetscReal nnmm;
  PetscFunctionBegin;
  info = MatNorm(pm,NORM_1,&nnmm);CHKERRQ(info);
  *nm=nnmm;
  PetscFunctionReturn(0);
}

#endif


