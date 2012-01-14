#include "tao_general.h"   /*I "tao_solver.h"  I*/
#include "taomat.h"
#include "taovec.h"


#undef __FUNCT__
#define __FUNCT__ "TaoMatDestroy"
/*@C
   TaoMatDestroy - Destroys the TaoMat object.

   Input Parameter:
.  TM - the matrix

   Level: beginner
@*/
int TaoMatDestroy( TaoMat* TM){
  TaoFunctionBegin;
  if (TM) {
    delete TM;
    TM=0;
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::Compatible"
/*@C
   Compatible - Confirms that the operation yy = this * xx is well defined.

   Input Parameters:
.  xx,yy -  vectors like those that will be used for matrix-vector multiplication

   Output Parameters:
.  flag - whether or not a matrix vector product can be performed.

   Level: developer
@*/
int TaoMat::Compatible(TaoVec *xx, TaoVec *yy, TaoTruth *flag){
  TaoFunctionBegin;
  *flag=TAO_FALSE;
  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoMat::GetDimensions"
/*@C
   GetDimensions - Gets the dimensions of the rowspace and columnspace of this matrix.

   Output Parameter:
+  m -  the number of rows
-  n - the number of columns

   Level: intermediate

@*/
int TaoMat::GetDimensions( TaoInt* m, TaoInt* n ){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::Multiply"
/*@C
   Multiply - Computes  ty = this * tx.

   Input Parameter:
.  tx -  the vector to be multiplied by this matrix.

   Output Parameter:
.  ty -  the destination vector


   Level: intermediate

.seealso TaoMat::MultiplyAdd(), TaoMat::MultiplyTranspose()

@*/
int TaoMat::Multiply(TaoVec* tx,TaoVec* ty){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}


#undef __FUNCT__
#define __FUNCT__ "TaoMat::MultiplyTranspose"
/*@C
   MultiplyTranspose - Multiplies the transpose of this matrix by a TaoVec.

   Input Parameter:
.  tx -  the vector to be multiplied by this matrix.

   Output Parameter:
.  ty -  the destination vector


   Level: intermediate

.seealso TaoMat::Multiply(), TaoMat::MultiplyTransposeAdd()
@*/
int TaoMat::MultiplyTranspose(TaoVec* tx,TaoVec* ty){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}


#undef __FUNCT__
#define __FUNCT__ "TaoMat::SetDiagonal"
/*@C
   SetDiagonal - Sets the diagonal elements of this matrix with the elements
   of the vector.

   Input Parameter:
.  tv -  the vector containing the diagonal elements

   Level: intermediate

.seealso TaoMat::AddDiagonal(),TaoMat::ShiftDiagonal()
@*/
int TaoMat::SetDiagonal(TaoVec* tv){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::AddDiagonal"
/*@C
   AddDiagonal - Adds the elements of the vector to the diagonal of this matrix.

   Input Parameter:
.  tv -  the vector containing the diagonal elements


   Level: intermediate

.seealso TaoMat::SetDiagonal(),TaoMat::ShiftDiagonal()
@*/
int TaoMat::AddDiagonal(TaoVec* tv){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::GetDiagonal"
/*@C
   GetDiagonal - Inserts the diagonal elements of this matrix into the vector.

   Output Parameter:
.  tv -  the vector containing the diagonal elements

   Level: intermediate

.seealso TaoMat::SetDiagonal()
@*/
int TaoMat::GetDiagonal(TaoVec* tv){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::ShiftDiagonal"
/*@C
   ShiftDiagonal - Adds this constant to the diagonal elements of this matrix.

   Input Parameter:
.  c -  the constant

   Level: intermediate

.seealso TaoMat::SetDiagonal(), TaoMat::AddDiagonal()
@*/
int TaoMat::ShiftDiagonal(double c){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::Presolve"
/*@C
   Presolve - Prepares to solve a linear system with this matrix by 
   doing some initial setup (e.g., computing parts of a preconditioner,
   such as matrix factorization).

   Input:

   Note:
   This routine is optional.  A linear solver object can also be used
   to implement this operation.

   Level: advanced

.seealso TaoMat::Solve()
@*/
int TaoMat::Presolve(){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::Solve"
/*@C
   Solve - Solves the linear system $this tx = tb$.

   Input Parameter:
.  tb -  right hand side

   Output Parameter:
+  tx -  solution vector
-  tt -  TAO_TRUE if solved to prescribed accuracy and TAO_FALSE otherwise

   Level: advanced

   Note:
   This routine is optional.  A linear solver object can also be used
   to implement this operation.

.seealso TaoApplication::GetLinearSolver(),TaoMat::Presolve()
@*/
int TaoMat::Solve(TaoVec* tb, TaoVec *tx, TaoTruth *tt){
  TaoFunctionBegin;
  SETERRQ(56,"No TAO linear solver has been set.");
  /* TaoFunctionReturn(1); */
}


#undef __FUNCT__
#define __FUNCT__ "TaoMat::CreateReducedMatrix"
/*@C
  CreateReducedMatrix - Creates a new matrix using the specified rows and columns of this matrix.

   Input Parameter:
+  row -  the rows of this matrix to be extracted
-  column -  the columns of this matrix to be extracted

   Output Parameter:
-  M - the new matrix


.seealso TaoMat::SetReducedMatrix()

   Level: intermediate
@*/
int TaoMat::CreateReducedMatrix(TaoIndexSet* row,TaoIndexSet* col, TaoMat **M){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
};

#undef __FUNCT__
#define __FUNCT__ "TaoMat::SetReducedMatrix"
/*@C
  SetReducedMatrix - Creates a new matrix using the specified rows and columns of this matrix.

   Input Parameter:
+  row -  the rows of this matrix to be extracted
-  col -  the columns of this matrix to be extracted

   Output Parameter:
.  M - the full new matrix


.seealso TaoMat::CreateReducedMatrix()

   Level: intermediate
@*/
int TaoMat::SetReducedMatrix(TaoMat *M, TaoIndexSet* row,TaoIndexSet* col){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
};


#undef __FUNCT__
#define __FUNCT__ "TaoMat::RowScale"
/*@C
  RowScale - Scales the rows of this matrix.

   Input Parameter:
.  scale - the scaling parameters

   Note:
   This method can also be called by using a pointer to TaoVec as
   input parameters.

   Level: intermediate
@*/
int TaoMat::RowScale(TaoVec* scale){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
};

#undef __FUNCT__
#define __FUNCT__ "TaoMat::ColScale"
/*@C
  ColScale - Scales the columns of this matrix.

   Input Parameter:
.  scale - the scaling parameters

   Level: intermediate
@*/
int TaoMat::ColScale(TaoVec* scale){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
};

#undef __FUNCT__
#define __FUNCT__ "fischer"
inline static double fischer(double a, double b)
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
#define __FUNCT__ "norm"
inline static double norm(double a, double b)
{
  return sqrt(a*a + b*b);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::D_Fischer"
/*@C
   D_Fischer - Calculates an element of the B-subdifferential of the 
   Fischer-Burmeister function for complementarity problems.
  
   Input Parameters:   
+  this - the jacobian of tf at tx
.  tx - current point
.  tf - function evaluated at tx
.  tl - lower bounds
.  tu - upper bounds
.  tt1 - work vector
-  tt2 - work vector

   Output Parameters:
+  tda - diagonal perturbation component of the result
-  tdb - row scaling component of the result

   Level: intermediate

.seealso TaoVec::Fischer()
@*/
int TaoMat::D_Fischer(TaoVec *tx, TaoVec *tf, TaoVec *tl, TaoVec *tu, 
		      TaoVec *tt1, TaoVec *tt2, TaoVec *tda, TaoVec *tdb )
{
  int i, info;
  TaoInt nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8;
  TaoTruth flag;
  TaoScalar *x,*f,*l,*u,*da,*db,*t1,*t2;
  TaoScalar ai,bi,ci,di,ei;

  TaoFunctionBegin;

  info = this->Compatible(tx, tx, &flag); CHKERRQ(info);
  if (flag==TAO_FALSE) {TaoFunctionReturn(1);}

  info = tx->GetArray(&x,&nn1);CHKERRQ(info);
  info = tf->GetArray(&f,&nn2);CHKERRQ(info);
  info = tl->GetArray(&l,&nn3);CHKERRQ(info);
  info = tu->GetArray(&u,&nn4);CHKERRQ(info);
  info = tda->GetArray(&da,&nn5);CHKERRQ(info);
  info = tdb->GetArray(&db,&nn6);CHKERRQ(info);
  info = tt1->GetArray(&t1,&nn7);CHKERRQ(info);
  info = tt2->GetArray(&t2,&nn8);CHKERRQ(info);

  if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4 || nn4!=nn5 || 
      nn5!=nn6 || nn6!=nn7 || nn7!=nn8) {
    TaoFunctionReturn(1);
  }

  for (i = 0; i < nn1; i++) {
    da[i] = 0;
    db[i] = 0;
    t1[i] = 0;

    if (TaoAbsScalar(f[i]) <= TAO_EPSILON) {
      if (l[i] > -TAO_INFINITY && TaoAbsScalar(x[i] - l[i]) <= TAO_EPSILON) {
        t1[i] = 1;
        da[i] = 1;
      }

      if (u[i] <  TAO_INFINITY && TaoAbsScalar(u[i] - x[i]) <= TAO_EPSILON) {
        t1[i] = 1;
        db[i] = 1;
      }
    }
  }

  info = tt1->RestoreArray(&t1,&nn7);CHKERRQ(info);
  info = tt2->RestoreArray(&t2,&nn8);CHKERRQ(info);
  info = this->Multiply(tt1, tt2); CHKERRQ(info);
  info = tt2->GetArray(&t2,&nn8); CHKERRQ(info);

  for (i = 0; i < nn1; i++) {
    if ((l[i] <= -TAO_INFINITY) && (u[i] >= TAO_INFINITY)) {
      da[i] = 0;
      db[i] = -1;
    } 
    else if (l[i] <= -TAO_INFINITY) {
      if (db[i] >= 1) {
        ai = norm(1, t2[i]);

        da[i] = -1/ai - 1;
        db[i] = -t2[i]/ai - 1;
      } 
      else {
        bi = u[i] - x[i];
        ai = norm(bi, f[i]);
        ai = TaoMax(TAO_EPSILON, ai);

        da[i] = bi / ai - 1;
        db[i] = -f[i] / ai - 1;
      }
    } 
    else if (u[i] >=  TAO_INFINITY) {
      if (da[i] >= 1) {
        ai = norm(1, t2[i]);

        da[i] = 1 / ai - 1;
        db[i] = t2[i] / ai - 1;
      } 
      else {
        bi = x[i] - l[i];
        ai = norm(bi, f[i]);
        ai = TaoMax(TAO_EPSILON, ai);

        da[i] = bi / ai - 1;
        db[i] = f[i] / ai - 1;
      }
    } 
    else if (l[i] == u[i]) {
      da[i] = -1;
      db[i] = 0;
    } 
    else {
      if (db[i] >= 1) {
        ai = norm(1, t2[i]);

        ci = 1 / ai + 1;
        di = t2[i] / ai + 1;
      } 
      else {
        bi = x[i] - u[i];
        ai = norm(bi, f[i]);
        ai = TaoMax(TAO_EPSILON, ai);

        ci = bi / ai + 1;
        di = f[i] / ai + 1;
      }

      if (da[i] >= 1) {
        bi = ci + di*t2[i];
        ai = norm(1, bi);

        bi = bi / ai - 1;
        ai = 1 / ai - 1;
      } 
      else {
        ei = fischer(u[i] - x[i], -f[i]);
        ai = norm(x[i] - l[i], ei);
        ai = TaoMax(TAO_EPSILON, ai);

        bi = ei / ai - 1;
        ai = (x[i] - l[i]) / ai - 1;
      }

      da[i] = ai + bi*ci;
      db[i] = bi*di;
    }
  }

  info = tda->RestoreArray(&da,&nn5);CHKERRQ(info);
  info = tdb->RestoreArray(&db,&nn6);CHKERRQ(info);

  info = tx->RestoreArray(&x,&nn1);CHKERRQ(info);
  info = tf->RestoreArray(&f,&nn2);CHKERRQ(info);
  info = tl->RestoreArray(&l,&nn3);CHKERRQ(info);
  info = tu->RestoreArray(&u,&nn4);CHKERRQ(info);
  info = tt2->RestoreArray(&t2,&nn8);CHKERRQ(info);
  TaoFunctionReturn(0);
};

#undef __FUNCT__
#define __FUNCT__ "sfischer"
inline static double sfischer(double a, double b, double c)
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
#define __FUNCT__ "snorm"
inline static double snorm(double a, double b, double c)
{
  return sqrt(a*a + b*b + 2.0*c*c);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::D_SFischer"
/*@C
   D_SFischer - Calculates an element of the B-subdifferential of the
   smoothed Fischer-Burmeister function for complementarity problems.
 
   Input Parameters: 
+  this - the jacobian of tf at tx
.  tx - current point
.  tf - function evaluated at tx
.  tl - lower bounds
.  tu - upper bounds
.  mu - smoothing parameter
.  tt1 - work vector
-  tt2 - work vector

   Output Parameter: 
+  tda - diagonal perturbation component of the result
.  tdb - row scaling component of the result
-  tdm - derivative with respect to scaling parameter

   Level: intermediate

.seealso TaoVec::SFischer()
@*/
int TaoMat::D_SFischer(TaoVec *tx, TaoVec *tf, 
                       TaoVec *tl, TaoVec *tu, double mu, 
                       TaoVec *tt1, TaoVec *tt2, 
                       TaoVec *tda, TaoVec *tdb, TaoVec *tdm)
{
  int i, info;
  TaoInt nn1, nn2, nn3, nn4, nn5, nn6, nn7;
  TaoTruth flag;
  TaoScalar *x, *f, *l, *u, *da, *db, *dm;
  TaoScalar ai, bi, ci, di, ei, fi;

  TaoFunctionBegin;

  if ((mu >= -TAO_EPSILON) && (mu <= TAO_EPSILON)) {
    tdm->SetToZero();
    D_Fischer(tx, tf, tl, tu, tt1, tt2, tda, tdb);
  } 
  else {
    info = this->Compatible(tx, tx, &flag); CHKERRQ(info);

    info = tx->GetArray(&x, &nn1); CHKERRQ(info);
    info = tf->GetArray(&f, &nn2); CHKERRQ(info);
    info = tl->GetArray(&l, &nn3); CHKERRQ(info);
    info = tu->GetArray(&u, &nn4); CHKERRQ(info);
    info = tda->GetArray(&da, &nn5); CHKERRQ(info);
    info = tdb->GetArray(&db, &nn6); CHKERRQ(info);
    info = tdm->GetArray(&dm, &nn7); CHKERRQ(info);

    if (nn1!=nn2 || nn2!=nn3 || nn3!=nn4 || nn4!=nn5 || nn5!=nn6 || nn6!=nn7) {
      TaoFunctionReturn(1);
    }

    for (i = 0; i < nn1; ++i) {
      if ((l[i] <= -TAO_INFINITY) && (u[i] >= TAO_INFINITY)) {
        da[i] = -mu;
        db[i] = -1;
        dm[i] = -x[i];
      } 
      else if (l[i] <= -TAO_INFINITY) {
        bi = u[i] - x[i];
        ai = snorm(bi, f[i], mu);
        ai = TaoMax(TAO_EPSILON, ai);

        da[i] = bi / ai - 1;
        db[i] = -f[i] / ai - 1;
        dm[i] = 2.0 * mu / ai;
      } 
      else if (u[i] >=  TAO_INFINITY) {
        bi = x[i] - l[i];
        ai = snorm(bi, f[i], mu);
        ai = TaoMax(TAO_EPSILON, ai);

        da[i] = bi / ai - 1;
        db[i] = f[i] / ai - 1;
        dm[i] = 2.0 * mu / ai;
      } 
      else if (l[i] == u[i]) {
        da[i] = -1;
        db[i] = 0;
        dm[i] = 0;
      } 
      else {
        bi = x[i] - u[i];
        ai = snorm(bi, f[i], mu);
        ai = TaoMax(TAO_EPSILON, ai);
  
        ci = bi / ai + 1;
        di = f[i] / ai + 1;
        fi = 2.0 * mu / ai;

        ei = sfischer(u[i] - x[i], -f[i], mu);
        ai = snorm(x[i] - l[i], ei, mu);
        ai = TaoMax(TAO_EPSILON, ai);
  
        bi = ei / ai - 1;
        ei = 2.0 * mu / ei;
        ai = (x[i] - l[i]) / ai - 1;
  
        da[i] = ai + bi*ci;
        db[i] = bi*di;
        dm[i] = ei + bi*fi;
      }
    }

    info = tx->RestoreArray(&x, &nn1); CHKERRQ(info);
    info = tf->RestoreArray(&f, &nn2); CHKERRQ(info);
    info = tl->RestoreArray(&l, &nn3); CHKERRQ(info);
    info = tu->RestoreArray(&u, &nn4); CHKERRQ(info);
    info = tda->RestoreArray(&da, &nn5); CHKERRQ(info);
    info = tdb->RestoreArray(&db, &nn6); CHKERRQ(info);
    info = tdm->RestoreArray(&dm, &nn7); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoMat::View"
/*@C
  View - Views the contents of this matrix.

   Input Parameter:

   Level: intermediate
@*/
int TaoMat::View(){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
};

#undef __FUNCT__
#define __FUNCT__ "TaoMat::Norm1"
/*@C
   Norm1 - Computes the 1-norm of the matrix.

   Output Parameter:
.  nm -  matrix 1-norm value

   Level: intermediate
@*/
int TaoMat::Norm1(double *nm){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

