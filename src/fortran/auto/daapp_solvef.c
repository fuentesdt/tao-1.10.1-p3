#include "petsc.h"
#include "petscfix.h"
/* daapp_solve.c */
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

#include "taodaapplication.h"
#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taodaappsolve_ PTAODAAPPSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taodaappsolve_ ptaodaappsolve
#else
#define taodaappsolve_ ptaodaappsolve_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taodaappsolve_ TAODAAPPSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taodaappsolve_ taodaappsolve
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taodaappsolve_(TAO_APPLICATION *daapplication,TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoDAAppSolve(*daapplication,*tao);
}
#if defined(__cplusplus)
}
#endif
