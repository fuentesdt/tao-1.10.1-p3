#include "petsc.h"
#include "petscfix.h"
/* tao_util.c */
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
#define taoconverged_maxits_ PTAOCONVERGED_MAXITS
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define taoconverged_maxits_ ptaoconverged_maxits__
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define taoconverged_maxits_ ptaoconverged_maxits
#else
#define taoconverged_maxits_ ptaoconverged_maxits_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoconverged_maxits_ TAOCONVERGED_MAXITS
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define taoconverged_maxits_ taoconverged_maxits__
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define taoconverged_maxits_ taoconverged_maxits
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoconverged_default_ PTAOCONVERGED_DEFAULT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define taoconverged_default_ ptaoconverged_default__
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define taoconverged_default_ ptaoconverged_default
#else
#define taoconverged_default_ ptaoconverged_default_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoconverged_default_ TAOCONVERGED_DEFAULT
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define taoconverged_default_ taoconverged_default__
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define taoconverged_default_ taoconverged_default
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdefaultmeritfunction_ PTAOSETDEFAULTMERITFUNCTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdefaultmeritfunction_ ptaosetdefaultmeritfunction
#else
#define taosetdefaultmeritfunction_ ptaosetdefaultmeritfunction_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdefaultmeritfunction_ TAOSETDEFAULTMERITFUNCTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdefaultmeritfunction_ taosetdefaultmeritfunction
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoconverged_maxits_(TAO_SOLVER *tao,void*dummy, int *__ierr ){
*__ierr = TaoConverged_MaxIts(*tao,dummy);
}
void PETSC_STDCALL taoconverged_default_(TAO_SOLVER *tao,void*dummy, int *__ierr ){
*__ierr = TaoConverged_Default(*tao,dummy);
}
void PETSC_STDCALL taosetdefaultmeritfunction_(TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoSetDefaultMeritFunction(*tao);
}
#if defined(__cplusplus)
}
#endif
