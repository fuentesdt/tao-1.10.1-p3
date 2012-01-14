#include "petsc.h"
#include "petscfix.h"
/* daapp_mgrid.c */
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
#define daappusemultigrid_ PDAAPPUSEMULTIGRID
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappusemultigrid_ pdaappusemultigrid
#else
#define daappusemultigrid_ pdaappusemultigrid_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappusemultigrid_ DAAPPUSEMULTIGRID
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappusemultigrid_ daappusemultigrid
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetupmultigrid_ PDAAPPSETUPMULTIGRID
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetupmultigrid_ pdaappsetupmultigrid
#else
#define daappsetupmultigrid_ pdaappsetupmultigrid_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetupmultigrid_ DAAPPSETUPMULTIGRID
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetupmultigrid_ daappsetupmultigrid
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL daappusemultigrid_(TAO_APPLICATION *daapplication,PetscInt *coarselevel, int *__ierr ){
*__ierr = DAAppUseMultigrid(*daapplication,*coarselevel);
}
void PETSC_STDCALL daappsetupmultigrid_(TAO_APPLICATION *daapplication,PetscInt *coarselevel, int *__ierr ){
*__ierr = DAAppSetupMultigrid(*daapplication,*coarselevel);
}
#if defined(__cplusplus)
}
#endif
