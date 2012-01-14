
#include "QpGenTaoLinsys.h"
#include "Ma27Solver.h"
#include "PetscVector.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"

QpGenTaoLinsys::QpGenTaoLinsys(  QpGen * factory,
				       QpGenData * data,
				       LinearAlgebraPackage * la,
				       SparseSymMatrix * Mat_in,
				       DoubleLinearSolver * solver_in )
  : QpGenLinsys( factory, data, la ), solver(solver_in)
{
  SpReferTo( Mat, Mat_in );
}


void QpGenTaoLinsys::putXDiagonal( OoqpVector& xdiag_in )
{
  PetscVector & xdiag = dynamic_cast<PetscVector&>(xdiag_in);
  int ierr;
  double * dXdiag;
  ierr = VecGetArray( xdiag.pv, &dXdiag ); assert( 0 == ierr );
  SimpleVectorHandle sXdiag( new SimpleVector( dXdiag, xdiag.length() ) );

  Mat->atPutDiagonal( 0, *sXdiag );

  ierr = VecRestoreArray( xdiag.pv, &dXdiag ); assert( 0 == ierr );
}


void QpGenTaoLinsys::putZDiagonal( OoqpVector& zdiag_in )
{
  PetscVector & zdiag = dynamic_cast<PetscVector&>(zdiag_in);
  int ierr;
  double * dZdiag;
  ierr = VecGetArray( zdiag.pv, &dZdiag ); assert( 0 == ierr );
  SimpleVectorHandle sZdiag( new SimpleVector( dZdiag, zdiag.length() ) );

  Mat->atPutDiagonal( nx + my, *sZdiag );

  ierr = VecRestoreArray( zdiag.pv, &dZdiag ); assert( 0 == ierr );
}


void QpGenTaoLinsys::solveCompressed( OoqpVector & rhs_in )
{
  PetscVector & rhs = dynamic_cast<PetscVector&>(rhs_in);
  int ierr;
  double * drhs;
  
  ierr = VecGetArray( rhs.pv, &drhs ); assert( 0 == ierr );
  SimpleVectorHandle srhs( new SimpleVector( drhs, rhs.length() ) );

  solver->solve( *srhs );

  ierr = VecRestoreArray( rhs.pv, &drhs ); assert( 0 == ierr );
}  


void QpGenTaoLinsys::factor(Data *prob, Variables *vars)
{
  this->QpGenLinsys::factor( prob, vars );
  solver->matrixChanged();
}


QpGenTaoLinsys::~QpGenTaoLinsys()
{
  delete solver;
}

