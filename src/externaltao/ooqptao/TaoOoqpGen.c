#include "TaoOoqpGen.h"


#undef PETSCVECTOR_H
#include "PetscVector.h"
#include "src/petsctao/matrix/taomat_petsc.h"
#include "QpGenTao.h"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "GondzioSolver.h"
#include "Status.h"
#include "OoqpPetscMonitor.h"

#define HAVE_GETRUSAGE
#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif /*HAVE_GETRUSAGE*/

int TaoDestroy_OOQPGEN(TAO_SOLVER tao, void *data);

static int TaoGetSizes_OOQP( TAO_SOLVER tao, void *data,
			     int * nx,   int * my, int * mz,
			     int * nnzQ, int * nnzA, int * nnzC );


//-------------------------------------------------------------------------
//wraps up the iteration data & checks for the convergence using TaoMonitor
extern "C"
int TaoStatusWrap(void *data_in)
{
  int info;

  StatusData *data = static_cast<StatusData*>(data_in);
  if (!data) {
    SETERRQ(1,"Error casting data_in to StatusData");
  }

  TAO_SOLVER solver = static_cast<TAO_SOLVER>(data->ctx);
  TAO_OOQPGEN *p_ooqp = (TAO_OOQPGEN*) solver->data;
  QpGenData *prob = (QpGenData *) data->data;
  QpGenVars *vars = (QpGenVars *) data->vars;
  QpGenResiduals *resids = (QpGenResiduals *) data->resids;
  Solver * ooqpsolver = (Solver *) data->solver;
  double dnorm = data->dnorm;

  double f, fnorm, cnorm;
  TaoTerminateReason reason;

  cnorm = resids->residualNorm()/dnorm;
  fnorm = sqrt(data->mu * vars->nComplementaryVariables );

  f = prob->objectiveValue( vars );

  info =  TaoMonitor(solver,p_ooqp->iterate++,f,fnorm,cnorm,1.0,&reason);
  CHKERRQ( info );
  int ooqpresult;
//   = ooqpsolver->defaultStatus( prob, vars,
//  					      resids, 
//  					      data->i, data->mu,
//  					      data->level );
  if( reason > 0 ) {
    ooqpresult = SUCCESSFUL_TERMINATION;
  } else if ( reason < 0 ) {
    if( reason == TAO_DIVERGED_MAXITS ) {
      ooqpresult = MAX_ITS_EXCEEDED;
    } else {
      ooqpresult = UNKNOWN;
    }
  } else { // reason == 0 
    ooqpresult = NOT_FINISHED;
  }
  
  return ooqpresult;
}

#undef __FUNCT__
#define __FUNCT__ "TaoSetUp_OOQPGEN"
int TaoSetUp_OOQPGEN(TAO_SOLVER tao, void *data)
{
  int info;
  double objective;
  TaoVec * xx;
  TaoFunctionBegin;
  TaoTruth tao_flag;

  TAO_OOQPGEN *p_ooqp = static_cast<TAO_OOQPGEN*>( data);
  if (!p_ooqp) {
    SETERRQ(1,"Could not cast data to TAO_OOQPGEN");
  }

  // Compute and save the Hessian matrix
  info = TaoGetSolution(tao,&xx); CHKERRQ( info );
  info = xx->Clone(&p_ooqp->xl);
  info = xx->Clone(&p_ooqp->xu);
  info = xx->Clone(&p_ooqp->g);

  if (tao->setupcalled) {
    if (!p_ooqp->gamma) {
      SETERRQ(1,"Trying to access null data member p_ooqp->gamma");
    }

    info = p_ooqp->gamma->Compatible(xx, &tao_flag); CHKERRQ(info);

    if (tao_flag==TAO_TRUE){
      TaoFunctionReturn(0);
    } else {
      info = TaoDestroy_OOQPGEN(tao, p_ooqp);CHKERRQ(info);
    }
  }

  info = TaoSetGradientVector(tao,p_ooqp->g);CHKERRQ(info);
  info = TaoSetVariableBounds(tao,p_ooqp->xl,p_ooqp->xu);CHKERRQ(info);
  info = TaoEvaluateVariableBounds(tao,p_ooqp->xl,p_ooqp->xu); CHKERRQ(info);
  xx->SetToZero();

  info = TaoGetHessian(tao,&p_ooqp->Q); CHKERRQ( info );
  info = TaoComputeHessian(tao, xx, p_ooqp->Q); CHKERRQ( info );

//-------------------------------up to here-------------------------
  //create an instance of QpGenTao (problem itself) & QpGenData (problem data)

  {
    int nx, my, mz;
    int nnzQ, nnzA, nnzC;
    TaoGetSizes_OOQP( tao, data, &nx, &my, &mz, &nnzQ, &nnzA, &nnzC );
    
    p_ooqp->qp = new QpGenTao( nx, my, mz, nnzQ, nnzA, nnzC );
  }
  // Change tao names to OOQP names so we don't go batty.
  TaoVec * c, *xlow, *xupp;
  TaoMat * A = tao->jacobian, *C = tao->CA;
  TaoVec * b = tao->vfunc;
  TaoVec * clow = tao->RXL, *cupp = tao->RXU;
  
  info = TaoGetGradient(tao,&c); CHKERRQ( info );
  info = TaoComputeFunctionGradient(tao,xx,&objective,c); CHKERRQ( info );
  info = TaoGetVariableBounds(tao,&xlow,&xupp); CHKERRQ( info );
  p_ooqp->prob =
    (QpGenData*) p_ooqp->qp->makeData( c, p_ooqp->Q,
					 xlow, xupp, A, b, C, clow, cupp );

  // Create some dual variables
  info = xx->Clone( &p_ooqp->gamma ); CHKERRQ( info );
  info = xx->Clone( &p_ooqp->phi ); CHKERRQ( info );
  info = TaoSetStepDirectionVector(tao,0);CHKERRQ(info);

  info = b->Clone( &p_ooqp->y ); CHKERRQ( info );
  info = clow->Clone( &p_ooqp->z ); CHKERRQ( info );

  // create an instance of QpGenVars (variables)

  p_ooqp->vars =
    (QpGenVars*) p_ooqp->qp->makeVariables( p_ooqp->prob, xx,
      					    p_ooqp->y, p_ooqp->z,
      					    p_ooqp->gamma, p_ooqp->phi );

  //create an instance of Residuals

  p_ooqp->resid = (QpGenResiduals*) p_ooqp->qp->makeResiduals( p_ooqp->prob );

  //create the Gondzio solver

  p_ooqp->ms = new GondzioSolver( p_ooqp->qp, p_ooqp->prob );

  //add a monitor, if necessary

  PetscTruth flg;
  info=PetscOptionsHasName(0,"-tao_monitor",&flg);CHKERRQ( info );
  if (flg==PETSC_TRUE){
    p_ooqp->ms->monitorSelf();
  }

  //add a convergence checking routine to the monitor (encapsulated TaoMonitor)
  CStatus * stat = new CStatus( TaoStatusWrap, tao );
  p_ooqp->ms->useStatus( stat );

  TaoFunctionReturn(0);
}
//-------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TaoDestroy_OOQPGEN"
int TaoDestroy_OOQPGEN(TAO_SOLVER tao, void *data)
{
  int info;

  TaoFunctionBegin;

  TAO_OOQPGEN *p_ooqp = (TAO_OOQPGEN*) data;

  if (tao->setupcalled)
    {
      //destroy all the internal variables created by TaoSetup

      info = TaoVecDestroy( p_ooqp->xl   ); CHKERRQ( info );
      info = TaoVecDestroy( p_ooqp->xu   ); CHKERRQ( info );
      info = TaoVecDestroy( p_ooqp->g   ); CHKERRQ( info );

      info = TaoVecDestroy( p_ooqp->y   ); CHKERRQ( info );
      info = TaoVecDestroy( p_ooqp->z   ); CHKERRQ( info );
      
      info = TaoVecDestroy( p_ooqp->gamma  ); CHKERRQ( info );
      info = TaoVecDestroy( p_ooqp->phi ); CHKERRQ( info );

      //destroy the instantiated variables, residuals, solver, problem & data      

      delete p_ooqp->ms;
      delete p_ooqp->vars;  
      delete p_ooqp->resid;
      delete p_ooqp->prob;
      delete p_ooqp->qp;
    }

  info = TaoSetGradientVector(tao,0);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,0);CHKERRQ(info);
  info = TaoSetVariableBounds(tao,0,0);CHKERRQ(info);

  
  TaoFunctionReturn(0);
}
//-------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TaoSolve_OOQPGEN"
static int TaoSolve_OOQPGEN(TAO_SOLVER tao,void *data)
{
  TaoTerminateReason reason;
  TaoTruth success;

  TaoFunctionBegin;

  TAO_OOQPGEN *p_ooqp = (TAO_OOQPGEN*) data;

  p_ooqp->iterate=1;
#ifdef HAVE_GETRUSAGE
  rusage before;
  getrusage( RUSAGE_SELF, &before );
#endif
  
  p_ooqp->ms->solve(p_ooqp->prob,p_ooqp->vars,p_ooqp->resid);
#ifdef HAVE_GETRUSAGE
    rusage  after;
    getrusage( RUSAGE_SELF, &after );
    printf( "QP solved in %f seconds.\n", 
	    (double) (after.ru_utime.tv_sec - before.ru_utime.tv_sec)
	    + (after.ru_utime.tv_usec - before.ru_utime.tv_usec) / 1000000.0 );
#endif

  TaoFunctionReturn(0);
}
//-------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TaoView_OOQPGEN"
static int TaoView_OOQPGEN(TAO_SOLVER tao,void *data)
{
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}
//-------------------------------------------------------------------------
EXTERN_C_BEGIN
#undef __FUNCT__
#define __FUNCT__ "TaoCreate_OOQPGEN"
int TaoCreate_OOQPGEN(TAO_SOLVER tao)
{
  int info;

  TaoFunctionBegin;
  
  TAO_OOQPGEN *p_ooqp;
  info = TaoNew(TAO_OOQPGEN,&p_ooqp); CHKERRQ(info);
  info = PetscLogObjectMemory(tao,sizeof(TAO_OOQPGEN)); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_OOQPGEN,(void*)p_ooqp); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_OOQPGEN,TaoDestroy_OOQPGEN); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,0); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_OOQPGEN); CHKERRQ(info);

  info = TaoSetMaximumIterates(tao,2000); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,4000); CHKERRQ(info);
  info = TaoSetTolerances(tao,1e-8,1e-8,1.0e-8,0.0); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END

static int TaoGetSizes_OOQP( TAO_SOLVER tao, void *data,
			     int * nx,   int * my, int * mz,
			     int * nnzQ, int * nnzA, int * nnzC )
{
  TAO_OOQPGEN *p_ooqp = (TAO_OOQPGEN*) data;
  int info;
  MatInfo m_info;

  {
    Mat QQ;
    info = TaoMatGetPetscMat(p_ooqp->Q,&QQ); CHKERRQ( info );
    MatGetInfo(QQ,MAT_LOCAL,&m_info);
    *nnzQ  = (int)m_info.nz_used;
    *nx    = (int)m_info.rows_global;
  }
  {
    Mat JJ;
    info = TaoMatGetPetscMat(tao->jacobian,&JJ); CHKERRQ( info );
    info = MatGetInfo(JJ,MAT_LOCAL,&m_info); CHKERRQ( info );
    *nnzA  = (int)m_info.nz_used;
    *my   = (int)m_info.rows_global;
  }
  {
    Mat CA;
    info = TaoMatGetPetscMat(tao->CA,&CA); CHKERRQ( info );
    info = MatGetInfo(CA,MAT_LOCAL,&m_info); CHKERRQ( info );
    *nnzC  = (int)m_info.nz_used;
    *mz  = (int)m_info.rows_global;
  }
  return 0;
}

#include "PetscSpSymMatrix.h"

extern "C"
int newTaoQpGen( TaoVec * x,
		 TaoVec * y,
		 TaoVec * z,
		 int nnzQs[], int nnzAs[], int nnzCs[],
		 TaoVec ** c, TaoMat ** Q,
		 TaoVec **  xlow, TaoVec **  xupp, 
		 TaoMat **  A,     TaoVec ** b,
		 TaoMat **  C,
		 TaoVec **  clow, TaoVec **  cupp )
{
  int ierr;

  int nx, my, mz;
  { 
    // Get the sizes
    Vec px, py, pz;
    ierr = TaoVecGetPetscVec( x, &px ); CHKERRQ( ierr );
    ierr = VecGetSize( px, &nx ); CHKERRQ( ierr );

    ierr = TaoVecGetPetscVec( y, &py ); CHKERRQ( ierr );
    ierr = VecGetSize( py, &my ); CHKERRQ( ierr );

    ierr = TaoVecGetPetscVec( z, &pz ); CHKERRQ( ierr );
    ierr = VecGetSize( pz, &mz ); CHKERRQ( ierr );
  }

  // Things the same shape as x
  ierr =  x->Clone( c );     CHKERRQ( ierr );
  ierr =  x->Clone( xlow );  CHKERRQ( ierr );
  ierr =  x->Clone( xupp );  CHKERRQ( ierr );
  // Things the same shape as y
  ierr =  y->Clone( b );     CHKERRQ( ierr );
  // Things the same shape as z
  ierr =  z->Clone( clow );     CHKERRQ( ierr );
  ierr =  z->Clone( cupp );     CHKERRQ( ierr );
  
  // Now for the Matrices
  {
    Mat pQ, pA, pC;
    TaoMatPetsc *tQ, *tA, *tC;

    ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, nx, nx,
			    PETSC_DEFAULT, nnzQs, &pQ ); CHKERRQ(ierr);
    ierr = MatSetFromOptions( pQ ); CHKERRQ(ierr);
    // We want the mat to be valid even if it is empty
    ierr = MatAssemblyBegin( pQ, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
    ierr = MatAssemblyEnd( pQ, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);

    ierr = TaoWrapPetscMat( pQ, &tQ); CHKERRQ(ierr);
    *Q = tQ;
    PetscSpSymMatrixHandle myQ( new PetscSpSymMatrix( pQ ) );
   myQ->abmaxnorm();

    ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, my, nx,
			    PETSC_DEFAULT, nnzAs, &pA ); CHKERRQ(ierr);
    // We want the mat to be valid even if it is empty
    ierr = MatAssemblyBegin( pA, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
    ierr = MatAssemblyEnd( pA, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
    ierr = TaoWrapPetscMat( pA, &tA); CHKERRQ(ierr);
    *A = tA;

    ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, mz, nx,
			    PETSC_DEFAULT, nnzCs, &pC ); CHKERRQ(ierr);
    // We want the mat to be valid even if it is empty
    ierr = MatAssemblyBegin( pC, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
    ierr = MatAssemblyEnd( pC, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
    ierr = TaoWrapPetscMat( pC, &tC); CHKERRQ(ierr);
    *C = tC;
  }

  return 0;
}

extern "C"
int freeTaoQpGen( TaoVec ** c, TaoMat ** Q,
		  TaoVec **  xlow,  TaoVec **  xupp,
		  TaoMat **  A,     TaoVec ** b,
		  TaoMat **  C, TaoVec **  clow, TaoVec **  cupp )
{
  int ierr;
  ierr = TaoVecDestroy( *c );    CHKERRQ( ierr ); *c = 0;
  ierr = TaoVecDestroy( *xlow ); CHKERRQ( ierr ); *xlow = 0;
  ierr = TaoVecDestroy( *xupp ); CHKERRQ( ierr ); *xupp = 0;

  ierr = TaoVecDestroy( *b );    CHKERRQ( ierr ); *b = 0;
  ierr = TaoVecDestroy( *clow ); CHKERRQ( ierr ); *clow = 0;
  ierr = TaoVecDestroy( *cupp ); CHKERRQ( ierr ); *cupp = 0;

  ierr = TaoMatDestroy( *Q ); CHKERRQ( ierr ); *Q = 0;
  ierr = TaoMatDestroy( *A ); CHKERRQ( ierr ); *A = 0;
  ierr = TaoMatDestroy( *C ); CHKERRQ( ierr ); *C = 0;

  return 0;
}

//-------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "TaoRegisterOoqp"
/*@C
   TaoRegisterOoqp - Registers the ooqp package.  Necessary to use OOQP solver

   Not Collective

   Level: intermediate

Concepts: TAO_SOLVER, register

.seealso:  TaoRegisterDestroy(), TaoRegisterAll()
@*/

extern "C"
int TaoRegisterOoqp(void)
{
  char path[] = "${TAO_DIR}/lib/${PETSC_ARCH}/libooqptao.so";
  int info;

  TaoFunctionBegin;

  info = TaoRegisterDynamic("tao_ooqpgen",path,"TaoCreate_OOQPGEN",TaoCreate_OOQPGEN); CHKERRQ(info);

  TaoFunctionReturn(0);

}








