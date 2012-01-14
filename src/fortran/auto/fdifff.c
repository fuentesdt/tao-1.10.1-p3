#include "petsc.h"
#include "petscfix.h"
/* fdiff.c */
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

#include "taoapp.h"
#include "src/tao_impl.h"
#include "src/petsctao/include/taopetsc.h"
#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetcoloring_ PTAOAPPSETCOLORING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetcoloring_ ptaoappsetcoloring
#else
#define taoappsetcoloring_ ptaoappsetcoloring_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetcoloring_ TAOAPPSETCOLORING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetcoloring_ taoappsetcoloring
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoappsetcoloring_(TAO_APPLICATION *taoapp,ISColoring *coloring, int *__ierr ){
*__ierr = TaoAppSetColoring(*taoapp,*coloring);
}
#if defined(__cplusplus)
}
#endif
