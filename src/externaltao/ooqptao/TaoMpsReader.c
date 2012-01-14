
#include "cMpsReader.h"
#include "cQpGenSparse.h"
#include "MpsReader.h"
#include "TaoMpsReader.h"

#include "PetscSpSymMatrix.h"
#include "PetscSpGenMatrix.h"
#include "PetscVector.h"

#include "petscmat.h"
#include "tao.h"
#include "TaoOoqpGen.h"
#include "QpGenTao.h"

#include <assert.h>

// These just pass through to keep the interface naming scheme constistent
int newTaoMpsReader( char filename[], void ** reader )
{
  int ierr;

  *reader = newCMpsReader( filename, &ierr );

  return ierr;
}

void TaoMpsReaderGetSizes( void * reader_, int * nx, int * ny, int * nz )
{
  cMpsReaderGetSizes( reader_, nx, ny, nz );
}

void TaoMpsReaderGetNNZ(  void * reader_in,
			  int * nnzQ, int * nnzA, int * nnzC )
{
  MpsReader * reader = (MpsReader *) reader_in;
  reader->numbersOfNonZeros( nnzQ, nnzA, nnzC );
}

int freeTaoMpsReader( void ** reader )
{
  int ierr = 0;

  try {
    freeCMpsReader( reader );
  } 
  catch ( ... ) {
    ierr = 1;
  }
  return ierr;
}
 

static int TaoOoqpInfBounds(PetscVector & d_petsc, PetscVector &id_petsc,
			     int uplo )
{
  double val = ( uplo > 0 ) ? TAO_INFINITY : TAO_NINFINITY;
  int n = d_petsc.getSize();
  assert( n == id_petsc.getSize() );

  int ierr;
  double * d, *id;
  ierr = VecGetArray( d_petsc.pv, &d ); CHKERRQ( ierr );
  ierr = VecGetArray( id_petsc.pv, &id ); CHKERRQ( ierr );
  for( int i = 0; i < n; i++ ) {
    if( id[i] == 0 ) d[i] = val;
  }
  ierr = VecRestoreArray( id_petsc.pv, &id ); CHKERRQ( ierr );
  ierr = VecRestoreArray( d_petsc.pv, &d ); CHKERRQ( ierr );
  return 0;
}

int  TaoMpsReaderReadQpGen( void * cReader, double * f,
				 TaoVec * c_tao, TaoMat * Q_tao,
				 TaoVec *  xlow_tao, TaoVec *  xupp_tao, 
				 TaoMat *  A_tao,     TaoVec * b_tao,
				 TaoMat *  C_tao,
				 TaoVec *  clow_tao,  TaoVec *  cupp_tao )
{
  MpsReader * reader = (MpsReader *) cReader;
  enum { upper = 1, lower = -1 };
  int ierr;
  TaoVec *ixupp_tao, *ixlow_tao, *icupp_tao, *iclow_tao;

  ierr = xlow_tao->Clone( &ixlow_tao ); CHKERRQ( ierr );
  ierr = xupp_tao->Clone( &ixupp_tao ); CHKERRQ( ierr );
    
  ierr = clow_tao->Clone( &iclow_tao ); CHKERRQ( ierr );
  ierr = cupp_tao->Clone( &icupp_tao ); CHKERRQ( ierr );
  {
    PetscVectorHandle c( TaoVecAsOoqpVec( c_tao ) );
    
        
    PetscVectorHandle xlow( TaoVecAsOoqpVec( xlow_tao ) );
    PetscVectorHandle ixlow( TaoVecAsOoqpVec( ixlow_tao ) );
    PetscVectorHandle xupp( TaoVecAsOoqpVec( xupp_tao ) );
    PetscVectorHandle ixupp( TaoVecAsOoqpVec( ixupp_tao ) );
    
    PetscVectorHandle b( TaoVecAsOoqpVec( b_tao ) );
    
    PetscVectorHandle clow( TaoVecAsOoqpVec( clow_tao ) );
    PetscVectorHandle iclow( TaoVecAsOoqpVec( iclow_tao ) );
    PetscVectorHandle cupp( TaoVecAsOoqpVec( cupp_tao ) );
    PetscVectorHandle icupp( TaoVecAsOoqpVec( icupp_tao ) );
    
    Mat Q_petsc, A_petsc, C_petsc;
    ierr = TaoMatGetPetscMat( Q_tao,    &Q_petsc );    CHKERRQ( ierr );
    ierr = TaoMatGetPetscMat( A_tao,    &A_petsc );    CHKERRQ( ierr );
    ierr = TaoMatGetPetscMat( C_tao,    &C_petsc );    CHKERRQ( ierr );
    
    PetscSpSymMatrixHandle Q( new PetscSpSymMatrix( Q_petsc ) );
    PetscSpGenMatrixHandle A( new PetscSpGenMatrix( A_petsc ) );
    PetscSpGenMatrixHandle C( new PetscSpGenMatrix( C_petsc ) );

  Q->abmaxnorm();
    reader->readQpGen( *c, *Q, *xlow, *ixlow, *xupp, *ixupp,
		       *A, *b,
		       *C, *clow, *iclow, *cupp, *icupp, ierr );

    *f = reader->objconst();

    if( ierr == 0 ) {
      ierr = TaoOoqpInfBounds( *xlow, *ixlow, lower ); CHKERRQ( ierr );
      ierr = TaoOoqpInfBounds( *xupp, *ixupp, upper ); CHKERRQ( ierr );
      
      ierr = TaoOoqpInfBounds( *clow, *iclow, lower ); CHKERRQ( ierr );
      ierr = TaoOoqpInfBounds( *cupp, *icupp, upper ); CHKERRQ( ierr );
    }
  }

  ierr = TaoVecDestroy( ixupp_tao ); CHKERRQ( ierr );
  ierr = TaoVecDestroy( ixlow_tao ); CHKERRQ( ierr );

  ierr = TaoVecDestroy( icupp_tao ); CHKERRQ( ierr );
  ierr = TaoVecDestroy( iclow_tao ); CHKERRQ( ierr );

  return ierr;
}

