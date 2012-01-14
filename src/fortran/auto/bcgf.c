#include "petsc.h"
#include "petscfix.h"
/* bcg.c */
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
#define taobcgsetrestarttol_ PTAOBCGSETRESTARTTOL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taobcgsetrestarttol_ ptaobcgsetrestarttol
#else
#define taobcgsetrestarttol_ ptaobcgsetrestarttol_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taobcgsetrestarttol_ TAOBCGSETRESTARTTOL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taobcgsetrestarttol_ taobcgsetrestarttol
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taobcgsetmethod_ PTAOBCGSETMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taobcgsetmethod_ ptaobcgsetmethod
#else
#define taobcgsetmethod_ ptaobcgsetmethod_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taobcgsetmethod_ TAOBCGSETMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taobcgsetmethod_ taobcgsetmethod
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taobcgsetrestarttol_(TAO_SOLVER *tao,double *eta, int *__ierr ){
*__ierr = TaoBCGSetRestartTol(*tao,*eta);
}
void PETSC_STDCALL taobcgsetmethod_(TAO_SOLVER *tao,TAO_CGTYPES *method, int *__ierr ){
*__ierr = TaoBCGSetMethod(*tao,*method);
}
#if defined(__cplusplus)
}
#endif
