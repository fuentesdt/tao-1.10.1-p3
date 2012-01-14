#include "petsc.h"
#include "petscfix.h"
/* line.c */
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
#define taogetcurrentsteplength_ PTAOGETCURRENTSTEPLENGTH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogetcurrentsteplength_ ptaogetcurrentsteplength
#else
#define taogetcurrentsteplength_ ptaogetcurrentsteplength_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogetcurrentsteplength_ TAOGETCURRENTSTEPLENGTH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogetcurrentsteplength_ taogetcurrentsteplength
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taogetcurrentsteplength_(TAO_SOLVER *tao,double *steplength,TaoTruth *newsearch, int *__ierr ){
*__ierr = TaoGetCurrentStepLength(*tao,steplength,
	(TaoTruth* )PetscToPointer((newsearch) ));
}
#if defined(__cplusplus)
}
#endif
