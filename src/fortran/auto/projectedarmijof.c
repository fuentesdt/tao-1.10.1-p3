#include "petsc.h"
#include "petscfix.h"
/* projectedarmijo.c */
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

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taocreateprojectedarmijolinesearch_ PTAOCREATEPROJECTEDARMIJOLINESEARCH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taocreateprojectedarmijolinesearch_ ptaocreateprojectedarmijolinesearch
#else
#define taocreateprojectedarmijolinesearch_ ptaocreateprojectedarmijolinesearch_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taocreateprojectedarmijolinesearch_ TAOCREATEPROJECTEDARMIJOLINESEARCH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taocreateprojectedarmijolinesearch_ taocreateprojectedarmijolinesearch
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taocreatendprojectedarmijolinesearch_ PTAOCREATENDPROJECTEDARMIJOLINESEARCH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taocreatendprojectedarmijolinesearch_ ptaocreatendprojectedarmijolinesearch
#else
#define taocreatendprojectedarmijolinesearch_ ptaocreatendprojectedarmijolinesearch_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taocreatendprojectedarmijolinesearch_ TAOCREATENDPROJECTEDARMIJOLINESEARCH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taocreatendprojectedarmijolinesearch_ taocreatendprojectedarmijolinesearch
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taocreateprojectedarmijolinesearch_(TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoCreateProjectedArmijoLineSearch(*tao);
}
void PETSC_STDCALL taocreatendprojectedarmijolinesearch_(TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoCreateNDProjectedArmijoLineSearch(*tao);
}
#if defined(__cplusplus)
}
#endif
