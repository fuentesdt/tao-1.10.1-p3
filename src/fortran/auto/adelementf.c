#include "petsc.h"
#include "petscfix.h"
/* adelement.c */
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
#define daappsetadelementfunctiongradient_ PDAAPPSETADELEMENTFUNCTIONGRADIENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetadelementfunctiongradient_ pdaappsetadelementfunctiongradient
#else
#define daappsetadelementfunctiongradient_ pdaappsetadelementfunctiongradient_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetadelementfunctiongradient_ DAAPPSETADELEMENTFUNCTIONGRADIENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetadelementfunctiongradient_ daappsetadelementfunctiongradient
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL daappsetadelementfunctiongradient_(TAO_APPLICATION *daapplication,
      int (*funcgrad)(int[2],DERIV_TYPE[4],DERIV_TYPE*,void*),
      int *flops,void*ctx, int *__ierr ){
*__ierr = DAAppSetADElementFunctionGradient(*daapplication,funcgrad,*flops,ctx);
}
#if defined(__cplusplus)
}
#endif
