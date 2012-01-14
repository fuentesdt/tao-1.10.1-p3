#include "petsc.h"
#include "petscfix.h"
/* petsckernel.c */
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

#include "tao_general.h"
#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoprintstatement_ PTAOPRINTSTATEMENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoprintstatement_ ptaoprintstatement
#else
#define taoprintstatement_ ptaoprintstatement_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoprintstatement_ TAOPRINTSTATEMENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoprintstatement_ taoprintstatement
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoprintstatement_(TAO_SOLVER *tao, char *statement, int *__ierr ){
*__ierr = TaoPrintStatement(*tao,statement);
}
#if defined(__cplusplus)
}
#endif
