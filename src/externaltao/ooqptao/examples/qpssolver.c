#include "TaoOoqpGen.h"
#include "taovec.h"
#include "src/petsctao/vector/taovec_petsc.h"
#include "TaoMpsReader.h"
#include <iostream.h>


static  char help[]=
"This is TaoOoqp";


typedef struct {
  Mat A;
  Vec B;
  double c;
} AppCtx;

int FormFunctionGradient(TAO_SOLVER, Vec, double *,Vec,void *);
static int FormHessian(TAO_SOLVER,Vec,Mat *, Mat *, MatStructure *, void *);


#undef __FUNCT__
#define __FUNCT__ "main"
int main( int argc, char *argv[] )
{
  PetscTruth flg;
  int info;
  char filename[51];
  AppCtx myapp;

  PetscInitialize( &argc, &argv,(char *)0,help );
  TaoInitialize( &argc, &argv,(char *)0,help );


  info = PetscOptionsGetString(0,"-qps",filename,50,&flg); CHKERRQ(info);
  if (flg==PETSC_FALSE){
    return 1;
  }



//-----------------------------------here we try to read QP-----------------------------------
  int nx, my, mz, nnzQ, nnzA, nnzC;;

  void * reader;
  info = newTaoMpsReader( filename, &reader );


  if( info != 0 ) {
    cout << "Couldn't read file " << filename << endl 
         << "For what it is worth, the error number is " << info << endl;
    return 1;
  }


  TaoMpsReaderGetSizes( reader, &nx, &my, &mz );

  // Create some template TaoVecs so we can clone them
  TaoVecPetsc *x, *y, *z;
  Vec px, py, pz;
  info = VecCreateSeq(PETSC_COMM_SELF, nx, &px ); CHKERRQ(info);
  info = TaoWrapPetscVec(px, &x); CHKERRQ(info);

  info = VecCreateSeq(PETSC_COMM_SELF, my, &py ); CHKERRQ(info);
  info = TaoWrapPetscVec(py, &y); CHKERRQ(info);

  info = VecCreateSeq(PETSC_COMM_SELF, mz, &pz ); CHKERRQ(info);
  info = TaoWrapPetscVec(pz, &z); CHKERRQ(info);
  
  double objconst;
  TaoVec *c, *xlow, *xupp, *b, *clow, *cupp;
  TaoMat *Q, *A, *C;

  {
    int *nnzQs = 0, *nnzAs = 0, *nnzCs = 0;
    nnzQs = new int[nx];
    if( my > 0 ) nnzAs = new int[my];
    if( mz > 0 ) nnzCs = new int[mz];
    
    TaoMpsReaderGetNNZ( reader, nnzQs, nnzAs, nnzCs );
    info = newTaoQpGen( x, y, z, nnzQs, nnzAs, nnzCs,
			&c, &Q, &xlow, &xupp,
			&A, &b, &C, &clow, &cupp ); CHKERRQ(info);
    delete [] nnzQs; delete [] nnzAs; delete [] nnzCs;
  } // end scope of nnzQs, nnzAs, nnzCs


  info = TaoMpsReaderReadQpGen( reader, &objconst,
				c, Q, xlow, xupp, 
				A,  b, C, clow, cupp ); CHKERRQ(info);

  info = freeTaoMpsReader( &reader ); CHKERRQ(info);
  
  TaoTerminateReason reason;
  TAO_SOLVER bqp;        
  TAO_APPLICATION qpsapp;        

  info = TaoPetscApplicationCreate(PETSC_COMM_SELF,&qpsapp); CHKERRQ(info);


  /* Will want to put this somewhere else eventually 
   info = TaoRegisterDynamic("tao_ooqpgen","${TAO_DIR}/lib/${PETSC_ARCH}/libooqptao.so","TaoCreate_OOQPGEN",TaoCreate_OOQPGEN);
  CHKERRQ(info);

  */
  /* Required for tao to recognize tao_ooqpgen option */
  info = TaoRegisterOoqp();

  info = TaoCreate(MPI_COMM_SELF,"tao_ooqpgen",&bqp); CHKERRQ(info);



  { // set up myapp, we need the petsc structures here
    Vec c_petsc;  Mat Q_petsc; 
    info = TaoVecGetPetscVec( c, &c_petsc ); CHKERRQ( info );
    info = TaoMatGetPetscMat( Q, &Q_petsc ); CHKERRQ( info );
    myapp.A = Q_petsc;
    myapp.B= c_petsc;
    myapp.c=objconst;
  }
  // Create a clone of c to give to Tao
  TaoVec * g; 
  info = c->Clone( &g ); CHKERRQ( info );
  {
    Vec g_petsc, x_petsc, xlow_petsc, xupp_petsc;
    Vec b_petsc, clow_petsc, cupp_petsc;
    Mat Q_petsc, A_petsc, C_petsc;
    info = TaoVecGetPetscVec( g, &g_petsc ); CHKERRQ( info );
    info = TaoVecGetPetscVec( x, &x_petsc ); CHKERRQ( info );
    info = TaoVecGetPetscVec( b, &b_petsc ); CHKERRQ( info );

    info = TaoVecGetPetscVec( xlow, &xlow_petsc ); CHKERRQ( info );
    info = TaoVecGetPetscVec( xupp, &xupp_petsc ); CHKERRQ( info );
    info = TaoVecGetPetscVec( clow, &clow_petsc ); CHKERRQ( info );
    info = TaoVecGetPetscVec( cupp, &cupp_petsc ); CHKERRQ( info );

    info = TaoMatGetPetscMat( Q, &Q_petsc ); CHKERRQ( info );
    info = TaoMatGetPetscMat( A, &A_petsc ); CHKERRQ( info );
    info = TaoMatGetPetscMat( C, &C_petsc ); CHKERRQ( info );
    info = TaoSetPetscFunctionGradient(qpsapp, x_petsc, g_petsc,
				       FormFunctionGradient, (void*)&myapp);
    CHKERRQ(info);
    info = TaoSetPetscHessian(qpsapp, Q_petsc, Q_petsc,
			      FormHessian,(void*) &myapp); CHKERRQ(info);
    info = TaoSetPetscJacobian(qpsapp, A_petsc,
			       TAO_NULL, (void*)&myapp); CHKERRQ(info);
    info = TaoSetPetscConstraintsFunction(qpsapp, b_petsc,
					  TAO_NULL, (void*)&myapp);
    CHKERRQ(info);
    info = TaoSetPetscVariableBounds(qpsapp,
				     xlow_petsc, xupp_petsc); CHKERRQ(info);
    info = TaoSetInequalityConstraints(qpsapp,
				       clow_petsc, C_petsc, cupp_petsc);
    CHKERRQ(info);
  } // end of the scope of the petsc variables
  info = TaoSetApplication(bqp,qpsapp); CHKERRQ(info);

  info = TaoSetFromOptions(bqp); CHKERRQ(info);
  info = TaoSolve(bqp); CHKERRQ(info);
  info = TaoGetTerminationReason(bqp,&reason); CHKERRQ(info);
  if (reason <= 0)
    PetscPrintf(MPI_COMM_WORLD,
		"Try a different TAO method, adjust some parameters, "
		"or check the function evaluation routines\n");

  double fff,gnorm,cnorm,xdiff;
  int iterate;
  info = TaoGetIterationData(bqp,&iterate, &fff, &gnorm,
			     &cnorm, &xdiff, &reason);
  printf(" Iterates: %d,    Optimal Solution:  %10.8e \n",
	 iterate,fff+objconst);

  info = TaoDestroy(bqp); CHKERRQ(info);
 
  info = TaoVecDestroy( g ); CHKERRQ(info);
  info = TaoVecDestroy( x ); CHKERRQ(info);
  info = TaoVecDestroy( y ); CHKERRQ(info);
  info = TaoVecDestroy( z ); CHKERRQ(info);

  info = freeTaoQpGen( &c, &Q, &xlow, &xupp,
		       &A, &b, &C, &clow, &cupp ); CHKERRQ(info);

  PetscFinalize();
  TaoFinalize();

  return 0;
}
int FormFunctionGradient(TAO_SOLVER tao, Vec X, double *fcn,Vec G,void *ptr){
  AppCtx* user=(AppCtx*)ptr;
  int info;
  PetscFunctionBegin;
  *fcn=user->c;
  info = VecCopy(user->B,G);
  PetscFunctionReturn(0);
}

int FormHessian(TAO_SOLVER tao,Vec X,Mat *H, Mat *Hpre, MatStructure *flg, void *ptr){
  AppCtx* user=(AppCtx*)ptr;
  PetscFunctionBegin;
  *H=user->A;
  *Hpre=user->A;
  PetscFunctionReturn(0);
}

