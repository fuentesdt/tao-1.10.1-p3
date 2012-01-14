/*$Id$*/

#include "private/fortranimpl.h"
#include "tao_solver.h"


#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogetterminationreason_    TAOGETTERMINATIONREASON
#define taocreate_                  TAOCREATE
#define taosetmethod_               TAOSETMETHOD
#define taogetsolution_             TAOGETSOLUTION
#define taogetgradient_             TAOGETGRADIENT
#define taogetvariablebounds_       TAOGETVARIABLEBOUNDS
#define taosetlinesearch_           TAOSETLINESEARCH
#define taogetiterationdata_        TAOGETSOLUTIONSTATUS
#define taogetsolutionstatus_       TAOGETSOLUTIONSTATUS
#define taogetlinearsolver_         TAOGETLINEARSOLVER
#define taosetoptionsprefix_        TAOSETOPTIONSPREFIX
#define taoappendoptionsprefix_     TAOAPPENDOPTIONSPREFIX
#define taogetoptionsprefix_        TAOGETOPTIONSPREFIX

#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define taogetterminationreason_    taogetterminationreason
#define taocreate_                  taocreate
#define taosetmethod_               taosetmethod
#define taogetsolution_             taogetsolution
#define taogetgradient_             taogetgradient
#define taogetvariablebounds_       taogetvariablebounds
#define taosetlinesearch_           taosetlinesearch
#define taogetiterationdata_        taogetsolutionstatus
#define taogetsolutionstatus_       taogetsolutionstatus
#define taogetlinearsolver_         taogetlinearsolver
#define taosetoptionsprefix_        taosetoptionsprefix
#define taoappendoptionsprefix_     taoappendoptionsprefix
#define taogetoptionsprefix_        taogetoptionsprefix

#endif

EXTERN_C_BEGIN

void PETSC_STDCALL taocreate_(MPI_Comm *comm, CHAR type PETSC_MIXED_LEN(len1),TAO_SOLVER *outtao,int *ierr PETSC_END_LEN(len1) PETSC_END_LEN(len2)){
  char *t;
  PetscTruth flg1;

  FIXCHAR(type,len1,t);
  *ierr = PetscStrncmp(t,"",len1-1,&flg1);

  if (flg1==PETSC_FALSE){
      *ierr = TaoCreate(MPI_Comm_f2c(*(MPI_Fint *)&*comm), t,outtao);
  } else if (flg1==PETSC_TRUE){
    *ierr = TaoCreate(MPI_Comm_f2c(*(MPI_Fint *)&*comm), 0,outtao);
  }
  FREECHAR(type,t);
}

void PETSC_STDCALL taogetterminationreason_(TAO_SOLVER *tao,TaoTerminateReason *r,int *info)
{
  *info = TaoGetTerminationReason(*tao,r);
}

void PETSC_STDCALL taosetmethod_(TAO_SOLVER *tao,CHAR type PETSC_MIXED_LEN(len),
                                int *ierr PETSC_END_LEN(len))
{
  char *t;

  FIXCHAR(type,len,t);
  *ierr = TaoSetMethod(*tao,t);
  FREECHAR(type,t);
}


static void (*f5)(TAO_SOLVER*,TaoVec**,TaoVec**,TaoVec**,TaoVec**,double*,double*,double*,TaoInt*,void*,int*);

EXTERN_C_END
static int ourtaolinesearch(TAO_SOLVER tao,TaoVec* x,TaoVec* g ,TaoVec* dx,TaoVec* w,double *f,double *f_full, double *step,TaoInt *flag,void *ctx)
{
  int info = 0;
  (*f5)(&tao,&x,&g,&dx,&w,f,f_full,step,flag,ctx,&info);CHKERRQ(info);
  return 0;
}
EXTERN_C_BEGIN

void PETSC_STDCALL taosetlinesearch_(TAO_SOLVER *tao,
				     void (*setup)(TAO_SOLVER,void*),
				     void (*options)(TAO_SOLVER,void*),
				     void (*func)(TAO_SOLVER*,TaoVec**,TaoVec**,TaoVec* *,TaoVec**, 
						  double*, double*, double*, TaoInt*, void*,int*),
				     void (*view)(TAO_SOLVER,void*),
				     void (*destroy)(TAO_SOLVER,void*),
				     void *ctx,int *info){
  f5 = func;
  *info = TaoSetLineSearch(*tao,0,0,ourtaolinesearch,0,0,ctx);
  /*  
   *info = TaoSetLineSearch(*tao,setup,options,ourtaolinesearch,view,destroy,ctx);
   */
}

					     

/* ------------------------------------------------------------------------- */


void PETSC_STDCALL taogetsolution_(TAO_SOLVER *tao,TaoVec **X,int *info ){
  *info = TaoGetSolution(*tao,X);
}


void PETSC_STDCALL taogetgradient_(TAO_SOLVER *tao,TaoVec **G,int *info ){
  *info = TaoGetSolution(*tao,G);
}



void PETSC_STDCALL taogetvariablebounds_(TAO_SOLVER *tao,TaoVec** XL,TaoVec** XU, int *info ){
  *info = TaoGetVariableBounds(*tao,XL,XU);
}



void PETSC_STDCALL taogetlinearsolver_(TAO_SOLVER *tao,TaoLinearSolver **S,int *info ){
  *info = TaoGetLinearSolver(*tao,S);
}


void PETSC_STDCALL taogetsolutionstatus_(TAO_SOLVER *tao, TaoInt *it, double *f, double *fnorm, double *cnorm, double *xdiff, TaoTerminateReason *reason,int*info){
  *info=TaoGetSolutionStatus(*tao,it,f,fnorm,cnorm,xdiff,reason);
}

void PETSC_STDCALL taosetoptionsprefix_(TAO_SOLVER *tao, CHAR prefix PETSC_MIXED_LEN(len), int *ierr PETSC_END_LEN(len))
{
  char *t;
  FIXCHAR(prefix,len,t);
  *ierr = TaoSetOptionsPrefix(*tao,t);
  FREECHAR(prefix,t);
}

void PETSC_STDCALL taoappendoptionsprefix_(TAO_SOLVER *tao, CHAR prefix PETSC_MIXED_LEN(len), int *ierr PETSC_END_LEN(len))
{
  char *t;
  FIXCHAR(prefix,len,t);
  *ierr = TaoAppendOptionsPrefix(*tao,t);
  FREECHAR(prefix,t);
}

void PETSC_STDCALL taogetoptionsprefix_(TAO_SOLVER *tao, CHAR prefix PETSC_MIXED_LEN(len), int *ierr PETSC_END_LEN(len))
{
  const char *tname;
  *ierr = TaoGetOptionsPrefix(*tao,&tname);
  *ierr = PetscStrncpy(prefix,tname,len);
}



#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogetconvergencehistory_   TAOGETCONVERGENCEHISTORY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define taogetconvergencehistory_   taogetconvergencehistory
#endif

void PETSC_STDCALL taogetconvergencehistory_(TAO_SOLVER *tao,
					     TaoInt *na, int *info) {
  TaoInt *cits;
  PetscScalar *ca;
  *info = TaoGetConvergenceHistory(*tao,&ca,&cits,na);

}

EXTERN_C_END




