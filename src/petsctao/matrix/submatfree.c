#include "submatfree.h"                /*I  "mat.h"  I*/

int ISCreateComplement(IS, Vec, IS *);
int VecISSetToConstant(IS, PetscScalar, Vec);


#undef __FUNCT__  
#define __FUNCT__ "MatCreateSubMatrixFree"
/*@C
  MatCreateSubMatrixFree - Creates a reduced matrix by masking a
  full matrix.

   Collective on matrix

   Input Parameters:
+  mat - matrix of arbitrary type
.  RowMask - the rows that will be zero
-  ColMask - the columns that will be zero

   Output Parameters:
.  J - New matrix

   Notes: 
   The user provides the input data and is responsible for destroying
   this data after matrix J has been destroyed.  
 
   Level: developer

.seealso: MatCreate()
@*/
int MatCreateSubMatrixFree(Mat mat,IS RowMask, IS ColMask, Mat *J)
{
  MPI_Comm     comm=((PetscObject)mat)->comm;
  MatSubMatFreeCtx ctx;
  int          info;
  PetscInt mloc,nloc,m,n;
  PetscFunctionBegin;

  info = PetscNew(_p_MatSubMatFreeCtx,&ctx);CHKERRQ(info);

  ctx->A=mat;
  //  ctx->Row=Row;
  //  ctx->Col=Col;

  info = MatGetSize(mat,&m,&n);CHKERRQ(info);
  info = MatGetLocalSize(mat,&mloc,&nloc);CHKERRQ(info);

  info = VecCreateMPI(comm,nloc,n,&ctx->VC);CHKERRQ(info);
  //  info = ISCreateComplement(Col, ctx->VC, &ctx->ColComplement);CHKERRQ(info);
  //  ctx->RowComplement=ctx->ColComplement;
  ctx->VR=ctx->VC;
  info =  PetscObjectReference((PetscObject)mat);CHKERRQ(info);

  info=ISCreateComplement(RowMask,ctx->VC,&ctx->RowComplement);CHKERRQ(info);
  info=ISCreateComplement(ColMask,ctx->VC,&ctx->ColComplement);CHKERRQ(info);
  /*
  info =  PetscObjectReference((PetscObject)ctx->RowComplement);CHKERRQ(info);
  info =  PetscObjectReference((PetscObject)ctx->ColComplement);CHKERRQ(info);
  */
  info = MatCreateShell(comm,mloc,nloc,m,n,ctx,J);CHKERRQ(info);

  //  info = MatShellSetOperation(*J,MATOP_GET_ROW,(void(*)())MatGetRow_SMF);CHKERRQ(info);
  //  info = MatShellSetOperation(*J,MATOP_RESTORE_ROW,(void(*)())MatRestoreRow_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_MULT,(void(*)())MatMult_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_DESTROY,(void(*)())MatDestroy_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_VIEW,(void(*)())MatView_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_MULT_TRANSPOSE,(void(*)())MatMultTranspose_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_DIAGONAL_SET,(void(*)())MatDiagonalSet_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_SHIFT,(void(*)())MatShift_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_EQUAL,(void(*)())MatEqual_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_SCALE,(void(*)())MatScale_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_TRANSPOSE,(void(*)())MatTranspose_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_GET_DIAGONAL,(void(*)())MatGetDiagonal_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_GET_SUBMATRICES,(void(*)())MatGetSubMatrices_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_NORM,(void(*)())MatNorm_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_DUPLICATE,(void(*)())MatDuplicate_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_GET_SUBMATRIX,(void(*)())MatGetSubMatrix_SMF);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_GET_ROW_MAX,(void(*)())MatDuplicate_SMF);CHKERRQ(info);

  info = PetscLogObjectParent(mat,*J); CHKERRQ(info);

  PetscFunctionReturn(0);  
}

#undef __FUNCT__  
#define __FUNCT__ "MatSMFResetRowColumn"
int MatSMFResetRowColumn(Mat mat,IS RowMask,IS ColMask){
  MatSubMatFreeCtx ctx;
  int           info;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = PetscObjectReference((PetscObject)RowMask);CHKERRQ(info);
  info = PetscObjectReference((PetscObject)ColMask);CHKERRQ(info);
  info = ISDestroy(ctx->RowComplement);CHKERRQ(info);
  info = ISDestroy(ctx->ColComplement);CHKERRQ(info);
  ctx->ColComplement=ColMask;
  ctx->RowComplement=RowMask;
  PetscFunctionReturn(0);  
}

#undef __FUNCT__  
#define __FUNCT__ "MatMult_SMF"
int MatMult_SMF(Mat mat,Vec a,Vec y)
{
  MatSubMatFreeCtx ctx;
  int           info;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = VecCopy(a,ctx->VR);CHKERRQ(info);
  info = VecISSetToConstant(ctx->ColComplement,0.0,ctx->VR);CHKERRQ(info);
  info = MatMult(ctx->A,ctx->VR,y);CHKERRQ(info);
  info = VecISSetToConstant(ctx->RowComplement,0.0,y);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatMultTranspose_SMF"
int MatMultTranspose_SMF(Mat mat,Vec a,Vec y)
{
  MatSubMatFreeCtx ctx;
  int info;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = VecCopy(a,ctx->VC);CHKERRQ(info);
  info = VecISSetToConstant(ctx->RowComplement,0.0,ctx->VC);CHKERRQ(info);
  info = MatMultTranspose(ctx->A,ctx->VC,y);CHKERRQ(info);
  info = VecISSetToConstant(ctx->ColComplement,0.0,y);CHKERRQ(info);
  PetscFunctionReturn(0);
} 

#undef __FUNCT__  
#define __FUNCT__ "MatDiagonalSet_SMF"
int MatDiagonalSet_SMF(Mat M, Vec D,InsertMode is)
{
  MatSubMatFreeCtx ctx;
  int           info;

  PetscFunctionBegin;
  info = MatShellGetContext(M,(void **)&ctx);CHKERRQ(info);
  info = MatDiagonalSet(ctx->A,D,is);
  PetscFunctionReturn(0);
} 

#undef __FUNCT__  
#define __FUNCT__ "MatDestroy_SMF"
int MatDestroy_SMF(Mat mat)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info=MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  //  info=ISDestroy(ctx->Row);CHKERRQ(info);
  //  info=ISDestroy(ctx->Col);CHKERRQ(info);
  info =MatDestroy(ctx->A);CHKERRQ(info);
  info=ISDestroy(ctx->RowComplement);CHKERRQ(info);
  info=ISDestroy(ctx->ColComplement);CHKERRQ(info);
  info=VecDestroy(ctx->VC);CHKERRQ(info);
  info = PetscFree(ctx); CHKERRQ(info);
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "MatView_SMF"
int MatView_SMF(Mat mat,PetscViewer viewer)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = MatView(ctx->A,viewer);CHKERRQ(info);
  //  info = ISView(ctx->Row,viewer);CHKERRQ(info);
  //  info = ISView(ctx->Col,viewer);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatShift_SMF"
int MatShift_SMF(Mat Y, PetscScalar a)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(Y,(void **)&ctx);CHKERRQ(info);
  info = MatShift(ctx->A,a);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatDuplicate_SMF"
int MatDuplicate_SMF(Mat mat,MatDuplicateOption op,Mat *M)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = MatCreateSubMatrixFree(ctx->A,ctx->RowComplement,ctx->ColComplement,M);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatEqual_SMF"
int MatEqual_SMF(Mat A,Mat B,PetscTruth *flg)
{
  int          info;
  MatSubMatFreeCtx  ctx1,ctx2;
  PetscTruth flg1,flg2,flg3;

  PetscFunctionBegin;
  info = MatShellGetContext(A,(void **)&ctx1);CHKERRQ(info);
  info = MatShellGetContext(B,(void **)&ctx2);CHKERRQ(info);
  info = ISEqual(ctx1->RowComplement,ctx2->RowComplement,&flg2);CHKERRQ(info);
  info = ISEqual(ctx1->ColComplement,ctx2->ColComplement,&flg3);CHKERRQ(info);
  if (flg2==PETSC_FALSE || flg3==PETSC_FALSE){
    *flg=PETSC_FALSE;
  } else {
    info = MatEqual(ctx1->A,ctx2->A,&flg1);CHKERRQ(info);
    if (flg1==PETSC_FALSE){ *flg=PETSC_FALSE;} 
    else { *flg=PETSC_TRUE;} 
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatScale_SMF"
int MatScale_SMF(Mat mat, PetscScalar a)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = MatScale(ctx->A,a);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatTranspose_SMF"
int MatTranspose_SMF(Mat mat,Mat *B)
{
  PetscFunctionBegin;
  PetscFunctionReturn(1);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetDiagonal_SMF"
int MatGetDiagonal_SMF(Mat mat,Vec v)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = MatGetDiagonal(ctx->A,v);CHKERRQ(info);
  //  info = VecISSetToConstant(ctx->RowComplement,0.0,v);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatDiagonalSet_SMF"
int MatGetRowMax_SMF(Mat M, Vec D)
{
  MatSubMatFreeCtx ctx;
  int           info;

  PetscFunctionBegin;
  info = MatShellGetContext(M,(void **)&ctx);CHKERRQ(info);
  info = MatGetRowMax(ctx->A,D,PETSC_NULL);
  //  info = VecISSetToConstant(ctx->RowComplement,0.0,D);CHKERRQ(info);
  PetscFunctionReturn(0);
} 


#undef __FUNCT__  
#define __FUNCT__ "MatGetSubMatrices_SMF"
int MatGetSubMatrices_SMF(Mat A,int n, IS *irow,IS *icol,MatReuse scall,Mat **B)
{
  int info,i;

  PetscFunctionBegin;
  if (scall == MAT_INITIAL_MATRIX) {
    info = PetscMalloc( (n+1)*sizeof(Mat),B );CHKERRQ(info);
  }

  for ( i=0; i<n; i++ ) {
    info = MatGetSubMatrix_SMF(A,irow[i],icol[i],scall,&(*B)[i]);CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetSubMatrix_SMF"
int MatGetSubMatrix_SMF(Mat mat,IS isrow,IS iscol,MatReuse cll,
			Mat *newmat)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  if (newmat){
    info=MatDestroy(*newmat);CHKERRQ(info);
  }
  info = MatCreateSubMatrixFree(ctx->A,isrow,iscol, newmat);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetRow_SMF"
int MatGetRow_SMF(Mat mat,PetscInt row,PetscInt *ncols,const PetscInt **cols,const PetscScalar **vals)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = MatGetRow(ctx->A,row,ncols,cols,vals);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatRestoreRow_SMF"
int MatRestoreRow_SMF(Mat mat,PetscInt row,PetscInt *ncols,const PetscInt **cols,const PetscScalar **vals)
{
  int          info;
  MatSubMatFreeCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = MatRestoreRow(ctx->A,row,ncols,cols,vals);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatGetColumnVector_SMF"
int MatGetColumnVector_SMF(Mat mat,Vec Y, PetscInt col)
{
  int    info;
  MatSubMatFreeCtx  ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = MatGetColumnVector(ctx->A,Y,col);CHKERRQ(info);
  //  info = VecISSetToConstant(ctx->RowComplement,0.0,Y);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatConvert_SMF"
int MatConvert_SMF(Mat mat,MatType newtype,Mat *NewMat)
{
  int info;
  int size;
  MatSubMatFreeCtx  ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  MPI_Comm_size(((PetscObject)mat)->comm,&size);
  PetscFunctionReturn(1);
}

#undef __FUNCT__  
#define __FUNCT__ "MatNorm_SMF"
int MatNorm_SMF(Mat mat,NormType type,PetscReal *norm)
{
  int info;
  MatSubMatFreeCtx  ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);

  if (type == NORM_FROBENIUS) {
    *norm = 1.0;
  } else if (type == NORM_1 || type == NORM_INFINITY) {
    *norm = 1.0;
  } else {
    SETERRQ(PETSC_ERR_SUP,"No two norm");
  }
  PetscFunctionReturn(0);
}
