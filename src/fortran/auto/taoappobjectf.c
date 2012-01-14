#include "petsc.h"
#include "petscfix.h"
/* taoappobject.c */
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
#define taodestroyapplication_ PTAODESTROYAPPLICATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taodestroyapplication_ ptaodestroyapplication
#else
#define taodestroyapplication_ ptaodestroyapplication_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taodestroyapplication_ TAODESTROYAPPLICATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taodestroyapplication_ taodestroyapplication
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taodestroyapplication_(TaoApplication *myapp, int *__ierr ){
*__ierr = TaoDestroyApplication(
	(TaoApplication* )PetscToPointer((myapp) ));
}
#if defined(__cplusplus)
}
#endif
