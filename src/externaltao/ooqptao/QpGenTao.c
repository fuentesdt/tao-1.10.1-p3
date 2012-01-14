
#include "QpGenTao.h"
#include "QpGenTaoLinsys.h"
#include "QpGenData.h"
#include "QpGenVars.h"

#include "PetscLinearAlgebraPackage.h"
#include "PetscSpSymMatrix.h"
#include "PetscSpGenMatrix.h"
#include "SparseSymMatrixHandle.h"
#include "PetscVector.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "PetscIterativeSolver.h"
#include "Ma27Solver.h"

#include "tao.h"

PetscVector * TaoVecAsOoqpVec( TaoVec * v )
{
  Vec pv;
  int info;
  info = TaoVecGetPetscVec( v, &pv ); assert( 0 == info );
  
  return new PetscVector( pv );
}

static PetscVector * TaoVecMaybeNilAsOoqpVec( TaoVec * v )
{
  if( !v ) {
    return new PetscVector( 0 );
  } else {
    return TaoVecAsOoqpVec( v );
  }
}

//-------------------------------------------------------------------------
//  creates an index vector based on the bound petsc vector provided this
//  function is needed because ooqp does not recognize
//  TAO_INFINITY/TAO_NINFINITY, instead, it uses an index vector on top of
//  the bound vector to identify if the actual bound exists or not
static PetscVector * TaoOoqpIndexBounds( TaoVec * tb )
{
  int m, n;
  int ierr;
  Vec b, ib;
  
  ierr = TaoVecGetPetscVec( tb, &b ); assert( 0 == ierr );
  ierr = VecGetLocalSize( b,  &m  ); assert( 0 == ierr );
  ierr = VecGetSize( b, &n ); ; assert( 0 == ierr );
  PetscVectorHandle hib( new PetscVector( m, n ) );
  {
    Vec ib = hib->pv;
    double *a, *ia;

    ierr = VecGetArray(   b,  &a  ); assert( 0 == ierr );
    ierr = VecGetArray(  ib,  &ia ); assert( 0 == ierr );
    
    int i;
    for( i = 0; i < m; i++ ) {
      if( TAO_INFINITY  != a[i] &&
	  TAO_NINFINITY != a[i] ) {
	ia[i] = 1.0;
      } else {
	ia[i] = 0.0; 
	a[i]  = 0.0;
      }
    }
    ierr = VecRestoreArray( b,   &a  ); assert( 0 == ierr );
    ierr = VecRestoreArray( ib, &ia ); assert( 0 == ierr );
  }
  return SpAsPointer( hib );
}
//-------------------------------------------------------------------------

QpGenTao::QpGenTao( int nx, int my, int mz,
				      int nnzQ_in, int nnzA_in, int nnzC_in ) :
  QpGen( nx, my, mz ), nnzQ( nnzQ_in ), nnzA( nnzA_in ), nnzC( nnzC_in )
{
  la = PetscLinearAlgebraPackage::soleInstance();
}

Data * QpGenTao::makeData( TaoVec * c_in,     TaoMat * Q_in,
			   TaoVec * xlow_in,  TaoVec * xupp_in,  
			   TaoMat * A_in,     TaoVec * b_in,
			   TaoMat * C_in,     
			   TaoVec * clow_in,    TaoVec * cupp_in )  
{
  int info;

  PetscVectorHandle c( TaoVecAsOoqpVec( c_in ));
 
  PetscVectorHandle  xlow( TaoVecAsOoqpVec( xlow_in ) );
  PetscVectorHandle  ixlow( TaoOoqpIndexBounds( xlow_in )  );
  PetscVectorHandle  xupp( TaoVecAsOoqpVec( xupp_in  ) );
  PetscVectorHandle  ixupp( TaoOoqpIndexBounds( xupp_in  ) );

  PetscVectorHandle b( TaoVecMaybeNilAsOoqpVec( b_in ) );

  PetscVectorHandle  clow( TaoVecMaybeNilAsOoqpVec( clow_in  ));
  PetscVectorHandle  iclow( TaoOoqpIndexBounds( clow_in  ));
  PetscVectorHandle  cupp( TaoVecMaybeNilAsOoqpVec( cupp_in ));
  PetscVectorHandle  icupp( TaoOoqpIndexBounds( cupp_in ));

  PetscSpSymMatrixHandle Q;
  PetscSpGenMatrixHandle A, C;

  if( Q_in ) {
    Mat pQ;
    info = TaoMatGetPetscMat(Q_in,&pQ); assert( 0 == info );
    Q = PetscSpSymMatrixHandle( new  PetscSpSymMatrix( pQ ) );
  } else {
    Q = PetscSpSymMatrixHandle( new  PetscSpSymMatrix( nx, 0 ) );
  }
  Q->abmaxnorm();
  if( A_in ) {
    Mat pA;
    info = TaoMatGetPetscMat(A_in,&pA); assert( 0 == info );
    A =  PetscSpGenMatrixHandle ( new PetscSpGenMatrix( pA ) );
  } else {
    A = PetscSpGenMatrixHandle ( new PetscSpGenMatrix( my, nx, 0 ) );
  }
  if( C_in ) {
    Mat pC;
    info = TaoMatGetPetscMat(C_in,&pC); assert( 0 == info );
    C = PetscSpGenMatrixHandle( new PetscSpGenMatrix( pC ));
  } else {
    C = PetscSpGenMatrixHandle( new PetscSpGenMatrix( mz, nx, 0 ));
  }

  QpGenData * 
    data = new QpGenData( PetscLinearAlgebraPackage::soleInstance(),
			  c, Q, xlow, ixlow, xupp, ixupp,
			  A, b,
			  C, clow, iclow, cupp, icupp );

  //data->print();
  return data;
}

//  Data * QpGenTao::makeData()
//  {
//    return new QpGenData( la, nx, my, mz, nnzQ, nnzA, nnzC );
//  }

LinearSystem * QpGenTao::makeLinsys( Data * prob_in )
{
  QpGenData * prob = (QpGenData *) prob_in;
  int n = nx + my + mz;
 
  SparseSymMatrixHandle Mat( new SparseSymMatrix( n,n + nnzQ
                                                      + nnzA + nnzC ) );
 
  SimpleVectorHandle v( new SimpleVector(n) );
  v->setToZero();
  Mat->setToDiagonal(*v);

  prob->putQIntoAt( *Mat, 0, 0 );
  prob->putAIntoAt( *Mat, nx, 0);
  prob->putCIntoAt( *Mat, nx + my, 0 );

  Ma27Solver * solver = new Ma27Solver( Mat );
  return new QpGenTaoLinsys( this, prob, la, Mat,
			  solver );
}

void QpGenTao::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			  OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  PetscVector & rhs  = dynamic_cast<PetscVector&>( rhs_in );
  PetscVector & rhs1 = dynamic_cast<PetscVector&>( rhs1_in );
  PetscVector & rhs2 = dynamic_cast<PetscVector&>( rhs2_in );
  PetscVector & rhs3 = dynamic_cast<PetscVector&>( rhs3_in );

  double *v, *v1, *v2, *v3;
  // Only sequential operations are supported
  int ierr;
  ierr = VecGetArray(rhs.pv, &v   ); assert( ierr == 0 );
  ierr = VecGetArray(rhs1.pv, &v1 ); assert( ierr == 0 );
  ierr = VecGetArray(rhs2.pv, &v2 ); assert( ierr == 0 );
  ierr = VecGetArray(rhs3.pv, &v3 ); assert( ierr == 0 );

  memcpy( &v[0],       v1, nx * sizeof( double ) );
  if( my > 0 ) memcpy( &v[nx],      v2, my * sizeof( double ) );
  if( mz > 0 ) memcpy( &v[nx + my], v3, mz * sizeof( double ) );

  ierr = VecRestoreArray(rhs.pv, &v   ); assert( ierr == 0 );
  ierr = VecRestoreArray(rhs1.pv, &v1 ); assert( ierr == 0 );
  ierr = VecRestoreArray(rhs2.pv, &v2 ); assert( ierr == 0 );
  ierr = VecRestoreArray(rhs3.pv, &v3 ); assert( ierr == 0 );
}

void QpGenTao::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				      OoqpVector& z_in, OoqpVector& vars_in )
{
  PetscVector & vars  = dynamic_cast<PetscVector&>( vars_in );
  PetscVector & x =     dynamic_cast<PetscVector&>( x_in );
  PetscVector & y =     dynamic_cast<PetscVector&>( y_in );
  PetscVector & z =     dynamic_cast<PetscVector&>( z_in );
  
  double *v, *v1, *v2, *v3;
  // Only sequential operations are supported
  int ierr;
  ierr = VecGetArray(vars.pv, &v   ); assert( ierr == 0 );
  ierr = VecGetArray(x.pv,    &v1 ); assert( ierr == 0 );
  ierr = VecGetArray(y.pv,    &v2 ); assert( ierr == 0 );
  ierr = VecGetArray(z.pv,    &v3 ); assert( ierr == 0 );

  memcpy( v1, &v[0],         nx * sizeof( double ) );
  if ( my > 0 ) memcpy( v2, &v[nx],        my * sizeof( double ) );
  if ( mz > 0 ) memcpy( v3, &v[nx + my],   mz * sizeof( double ) );
  
  ierr = VecRestoreArray(vars.pv,  &v   ); assert( ierr == 0 );
  ierr = VecRestoreArray(x.pv,     &v1 ); assert( ierr == 0 );
  ierr = VecRestoreArray(y.pv,     &v2 ); assert( ierr == 0 );
  ierr = VecRestoreArray(z.pv,     &v3 ); assert( ierr == 0 );

}

Variables     * QpGenTao:: makeVariables( Data * problem, 
                                          TaoVec * x_in,
					  TaoVec * y_in,
					  TaoVec * z_in,
					  TaoVec * gamma_in,
					  TaoVec * phi_in )
{
  int info;
  QpGenData * prob = (QpGenData *) problem;

  PetscVectorHandle x( TaoVecAsOoqpVec( x_in ) );
  PetscVectorHandle y( TaoVecAsOoqpVec( y_in ) );
  PetscVectorHandle z( TaoVecAsOoqpVec( z_in ) );
  PetscVectorHandle gamma( TaoVecAsOoqpVec( gamma_in ) );
  PetscVectorHandle phi( TaoVecAsOoqpVec( phi_in ) );

  PetscVectorHandle v( new PetscVector( nx ) );
  PetscVectorHandle w( new PetscVector( nx ) );

  PetscVectorHandle s( new PetscVector( mz ) );
  PetscVectorHandle t( new PetscVector( mz ) );
  PetscVectorHandle lambda( new PetscVector( mz ) );
  PetscVectorHandle u( new PetscVector( mz ) );
  PetscVectorHandle pi( new PetscVector( mz ) );

  return new QpGenVars( x, s, y, z, v,
			gamma, w, phi, t, lambda,
			u, pi, 
			prob->ixlow, prob->ixupp,
			prob->iclow, prob->icupp );
}


