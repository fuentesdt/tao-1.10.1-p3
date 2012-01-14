#include "petsc.h"
#include "petscfix.h"
/* tao_app_util.c */
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
#define taoappsetvariablebounds_ PTAOAPPSETVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetvariablebounds_ ptaoappsetvariablebounds
#else
#define taoappsetvariablebounds_ ptaoappsetvariablebounds_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetvariablebounds_ TAOAPPSETVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetvariablebounds_ taoappsetvariablebounds
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoappsetvariablebounds_(TAO_APPLICATION *taoapp,Vec XL,Vec XU, int *__ierr ){
*__ierr = TaoAppSetVariableBounds(*taoapp,
	(Vec)PetscToPointer((XL) ),
	(Vec)PetscToPointer((XU) ));
}
#if defined(__cplusplus)
}
#endif
