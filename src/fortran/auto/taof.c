#include "petsc.h"
#include "petscfix.h"
/* tao.c */
/* Fortran interface file */

/*
* This file was generated automatically by bfort from the C source
* file.  
 */

#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" { 
#endif 
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
} 
#endif 

#else

#define PetscToPointer(a) (*(long *)(a))
#define PetscFromPointer(a) (long)(a)
#define PetscRmPointer(a)
#endif

#include "tao_solver.h"
#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoview_ PTAOVIEW
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoview_ ptaoview
#else
#define taoview_ ptaoview_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoview_ TAOVIEW
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoview_ taoview
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetgradienttolerances_ PTAOSETGRADIENTTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetgradienttolerances_ ptaosetgradienttolerances
#else
#define taosetgradienttolerances_ ptaosetgradienttolerances_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetgradienttolerances_ TAOSETGRADIENTTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetgradienttolerances_ taosetgradienttolerances
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogetgradienttolerances_ PTAOGETGRADIENTTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogetgradienttolerances_ ptaogetgradienttolerances
#else
#define taogetgradienttolerances_ ptaogetgradienttolerances_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogetgradienttolerances_ TAOGETGRADIENTTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogetgradienttolerances_ taogetgradienttolerances
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetfunctionlowerbound_ PTAOSETFUNCTIONLOWERBOUND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetfunctionlowerbound_ ptaosetfunctionlowerbound
#else
#define taosetfunctionlowerbound_ ptaosetfunctionlowerbound_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetfunctionlowerbound_ TAOSETFUNCTIONLOWERBOUND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetfunctionlowerbound_ taosetfunctionlowerbound
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetmaximumiterates_ PTAOSETMAXIMUMITERATES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetmaximumiterates_ ptaosetmaximumiterates
#else
#define taosetmaximumiterates_ ptaosetmaximumiterates_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetmaximumiterates_ TAOSETMAXIMUMITERATES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetmaximumiterates_ taosetmaximumiterates
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetmaximumfunctionevaluations_ PTAOSETMAXIMUMFUNCTIONEVALUATIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetmaximumfunctionevaluations_ ptaosetmaximumfunctionevaluations
#else
#define taosetmaximumfunctionevaluations_ ptaosetmaximumfunctionevaluations_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetmaximumfunctionevaluations_ TAOSETMAXIMUMFUNCTIONEVALUATIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetmaximumfunctionevaluations_ taosetmaximumfunctionevaluations
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosettolerances_ PTAOSETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosettolerances_ ptaosettolerances
#else
#define taosettolerances_ ptaosettolerances_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosettolerances_ TAOSETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosettolerances_ taosettolerances
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogettolerances_ PTAOGETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogettolerances_ ptaogettolerances
#else
#define taogettolerances_ ptaogettolerances_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogettolerances_ TAOGETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogettolerances_ taogettolerances
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoclearmonitor_ PTAOCLEARMONITOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoclearmonitor_ ptaoclearmonitor
#else
#define taoclearmonitor_ ptaoclearmonitor_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoclearmonitor_ TAOCLEARMONITOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoclearmonitor_ taoclearmonitor
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetconvergencehistory_ PTAOSETCONVERGENCEHISTORY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetconvergencehistory_ ptaosetconvergencehistory
#else
#define taosetconvergencehistory_ ptaosetconvergencehistory_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetconvergencehistory_ TAOSETCONVERGENCEHISTORY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetconvergencehistory_ taosetconvergencehistory
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosolve_ PTAOSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosolve_ ptaosolve
#else
#define taosolve_ ptaosolve_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosolve_ TAOSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosolve_ taosolve
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosettrustregiontolerance_ PTAOSETTRUSTREGIONTOLERANCE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosettrustregiontolerance_ ptaosettrustregiontolerance
#else
#define taosettrustregiontolerance_ ptaosettrustregiontolerance_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosettrustregiontolerance_ TAOSETTRUSTREGIONTOLERANCE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosettrustregiontolerance_ taosettrustregiontolerance
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogettrustregionradius_ PTAOGETTRUSTREGIONRADIUS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogettrustregionradius_ ptaogettrustregionradius
#else
#define taogettrustregionradius_ ptaogettrustregionradius_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogettrustregionradius_ TAOGETTRUSTREGIONRADIUS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogettrustregionradius_ taogettrustregionradius
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosettrustregionradius_ PTAOSETTRUSTREGIONRADIUS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosettrustregionradius_ ptaosettrustregionradius
#else
#define taosettrustregionradius_ ptaosettrustregionradius_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosettrustregionradius_ TAOSETTRUSTREGIONRADIUS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosettrustregionradius_ taosettrustregionradius
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetvariablebounds_ PTAOSETVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetvariablebounds_ ptaosetvariablebounds
#else
#define taosetvariablebounds_ ptaosetvariablebounds_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetvariablebounds_ TAOSETVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetvariablebounds_ taosetvariablebounds
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogetdualvariables_ PTAOGETDUALVARIABLES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogetdualvariables_ ptaogetdualvariables
#else
#define taogetdualvariables_ ptaogetdualvariables_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogetdualvariables_ TAOGETDUALVARIABLES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogetdualvariables_ taogetdualvariables
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetterminationreason_ PTAOSETTERMINATIONREASON
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetterminationreason_ ptaosetterminationreason
#else
#define taosetterminationreason_ ptaosetterminationreason_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetterminationreason_ TAOSETTERMINATIONREASON
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetterminationreason_ taosetterminationreason
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoview_(TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoView(*tao);
}
void PETSC_STDCALL taosetgradienttolerances_(TAO_SOLVER *tao,double *gatol,double *grtol,double *gttol, int *__ierr ){
*__ierr = TaoSetGradientTolerances(*tao,*gatol,*grtol,*gttol);
}
void PETSC_STDCALL taogetgradienttolerances_(TAO_SOLVER *tao,double *gatol,double *grtol,double *gttol, int *__ierr ){
*__ierr = TaoGetGradientTolerances(*tao,gatol,grtol,gttol);
}
void PETSC_STDCALL taosetfunctionlowerbound_(TAO_SOLVER *tao,double *fmin, int *__ierr ){
*__ierr = TaoSetFunctionLowerBound(*tao,*fmin);
}
void PETSC_STDCALL taosetmaximumiterates_(TAO_SOLVER *tao,TaoInt *maxits, int *__ierr ){
*__ierr = TaoSetMaximumIterates(*tao,*maxits);
}
void PETSC_STDCALL taosetmaximumfunctionevaluations_(TAO_SOLVER *tao,TaoInt *nfcn, int *__ierr ){
*__ierr = TaoSetMaximumFunctionEvaluations(*tao,*nfcn);
}
void PETSC_STDCALL taosettolerances_(TAO_SOLVER *tao,double *fatol,double *frtol,double *catol,double *crtol, int *__ierr ){
*__ierr = TaoSetTolerances(*tao,*fatol,*frtol,*catol,*crtol);
}
void PETSC_STDCALL taogettolerances_(TAO_SOLVER *tao,double *fatol,double *frtol,double *catol,double *crtol, int *__ierr ){
*__ierr = TaoGetTolerances(*tao,fatol,frtol,catol,crtol);
}
void PETSC_STDCALL taoclearmonitor_(TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoClearMonitor(*tao);
}
void PETSC_STDCALL taosetconvergencehistory_(TAO_SOLVER *tao,double *a,TaoInt *its,TaoInt *na,TaoTruth *reset, int *__ierr ){
*__ierr = TaoSetConvergenceHistory(*tao,a,
	(TaoInt* )PetscToPointer((its) ),*na,*reset);
}
void PETSC_STDCALL taosolve_(TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoSolve(*tao);
}
void PETSC_STDCALL taosettrustregiontolerance_(TAO_SOLVER *tao,double *steptol, int *__ierr ){
*__ierr = TaoSetTrustRegionTolerance(*tao,*steptol);
}
void PETSC_STDCALL taogettrustregionradius_(TAO_SOLVER *tao,double *radius, int *__ierr ){
*__ierr = TaoGetTrustRegionRadius(*tao,radius);
}
void PETSC_STDCALL taosettrustregionradius_(TAO_SOLVER *tao,double *radius, int *__ierr ){
*__ierr = TaoSetTrustRegionRadius(*tao,*radius);
}
void PETSC_STDCALL taosetvariablebounds_(TAO_SOLVER *tao,TaoVec *xxll,TaoVec *xxuu, int *__ierr ){
*__ierr = TaoSetVariableBounds(*tao,
	(TaoVec* )PetscToPointer((xxll) ),
	(TaoVec* )PetscToPointer((xxuu) ));
}
void PETSC_STDCALL taogetdualvariables_(TAO_SOLVER *tao,TaoVec *DXL,TaoVec *DXU, int *__ierr ){
*__ierr = TaoGetDualVariables(*tao,
	(TaoVec* )PetscToPointer((DXL) ),
	(TaoVec* )PetscToPointer((DXU) ));
}
void PETSC_STDCALL taosetterminationreason_(TAO_SOLVER *tao,TaoTerminateReason *reason, int *__ierr ){
*__ierr = TaoSetTerminationReason(*tao,*reason);
}
#if defined(__cplusplus)
}
#endif
