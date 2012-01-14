#include "petsc.h"
#include "petscfix.h"
/* tao_mfjdef.c */
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
#define mattaomfdefaultsetumin_ PMATTAOMFDEFAULTSETUMIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfdefaultsetumin_ pmattaomfdefaultsetumin
#else
#define mattaomfdefaultsetumin_ pmattaomfdefaultsetumin_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfdefaultsetumin_ MATTAOMFDEFAULTSETUMIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfdefaultsetumin_ mattaomfdefaultsetumin
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL mattaomfdefaultsetumin_(Mat A,double *umin, int *__ierr ){
*__ierr = MatTaoMFDefaultSetUmin(
	(Mat)PetscToPointer((A) ),*umin);
}
#if defined(__cplusplus)
}
#endif
