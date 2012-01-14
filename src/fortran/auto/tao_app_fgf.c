#include "petsc.h"
#include "petscfix.h"
/* tao_app_fg.c */
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

#include "tao.h"
#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappcomputeobjective_ PTAOAPPCOMPUTEOBJECTIVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappcomputeobjective_ ptaoappcomputeobjective
#else
#define taoappcomputeobjective_ ptaoappcomputeobjective_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappcomputeobjective_ TAOAPPCOMPUTEOBJECTIVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappcomputeobjective_ taoappcomputeobjective
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappcomputegradient_ PTAOAPPCOMPUTEGRADIENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappcomputegradient_ ptaoappcomputegradient
#else
#define taoappcomputegradient_ ptaoappcomputegradient_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappcomputegradient_ TAOAPPCOMPUTEGRADIENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappcomputegradient_ taoappcomputegradient
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappcomputeobjectiveandgradient_ PTAOAPPCOMPUTEOBJECTIVEANDGRADIENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappcomputeobjectiveandgradient_ ptaoappcomputeobjectiveandgradient
#else
#define taoappcomputeobjectiveandgradient_ ptaoappcomputeobjectiveandgradient_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappcomputeobjectiveandgradient_ TAOAPPCOMPUTEOBJECTIVEANDGRADIENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappcomputeobjectiveandgradient_ taoappcomputeobjectiveandgradient
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsethessianmat_ PTAOAPPSETHESSIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsethessianmat_ ptaoappsethessianmat
#else
#define taoappsethessianmat_ ptaoappsethessianmat_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsethessianmat_ TAOAPPSETHESSIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsethessianmat_ taoappsethessianmat
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappgethessianmat_ PTAOAPPGETHESSIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappgethessianmat_ ptaoappgethessianmat
#else
#define taoappgethessianmat_ ptaoappgethessianmat_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappgethessianmat_ TAOAPPGETHESSIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappgethessianmat_ taoappgethessianmat
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoapphessiansolve_ PTAOAPPHESSIANSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoapphessiansolve_ ptaoapphessiansolve
#else
#define taoapphessiansolve_ ptaoapphessiansolve_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoapphessiansolve_ TAOAPPHESSIANSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoapphessiansolve_ taoapphessiansolve
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappresetcounters_ PTAOAPPRESETCOUNTERS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappresetcounters_ ptaoappresetcounters
#else
#define taoappresetcounters_ ptaoappresetcounters_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappresetcounters_ TAOAPPRESETCOUNTERS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappresetcounters_ taoappresetcounters
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappcounters_ PTAOAPPCOUNTERS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappcounters_ ptaoappcounters
#else
#define taoappcounters_ ptaoappcounters_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappcounters_ TAOAPPCOUNTERS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappcounters_ taoappcounters
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoappcomputeobjective_(TAO_APPLICATION *taoapp,Vec X,double *f, int *__ierr ){
*__ierr = TaoAppComputeObjective(*taoapp,
	(Vec)PetscToPointer((X) ),f);
}
void PETSC_STDCALL taoappcomputegradient_(TAO_APPLICATION *taoapp,Vec X,Vec G, int *__ierr ){
*__ierr = TaoAppComputeGradient(*taoapp,
	(Vec)PetscToPointer((X) ),
	(Vec)PetscToPointer((G) ));
}
void PETSC_STDCALL taoappcomputeobjectiveandgradient_(TAO_APPLICATION *taoapp,Vec X,double *f,Vec G, int *__ierr ){
*__ierr = TaoAppComputeObjectiveAndGradient(*taoapp,
	(Vec)PetscToPointer((X) ),f,
	(Vec)PetscToPointer((G) ));
}
void PETSC_STDCALL taoappsethessianmat_(TAO_APPLICATION *taoapp,Mat H,Mat HP, int *__ierr ){
*__ierr = TaoAppSetHessianMat(*taoapp,
	(Mat)PetscToPointer((H) ),
	(Mat)PetscToPointer((HP) ));
}
void PETSC_STDCALL taoappgethessianmat_(TAO_APPLICATION *taoapp,Mat *H,Mat *HP, int *__ierr ){
*__ierr = TaoAppGetHessianMat(*taoapp,H,HP);
}
void PETSC_STDCALL taoapphessiansolve_(TAO_APPLICATION *taoapp,Vec Vin,Vec Vout,PetscTruth *success, int *__ierr ){
*__ierr = TaoAppHessianSolve(*taoapp,
	(Vec)PetscToPointer((Vin) ),
	(Vec)PetscToPointer((Vout) ),success);
}
void PETSC_STDCALL taoappresetcounters_(TAO_APPLICATION *taoapp, int *__ierr ){
*__ierr = TaoAppResetCounters(*taoapp);
}
void PETSC_STDCALL taoappcounters_(TAO_APPLICATION *taoapp,TaoInt stats[4], int *__ierr ){
*__ierr = TaoAppCounters(*taoapp,
	(TaoInt* )PetscToPointer((stats) ));
}
#if defined(__cplusplus)
}
#endif
