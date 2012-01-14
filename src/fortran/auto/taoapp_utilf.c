#include "petsc.h"
#include "petscfix.h"
/* taoapp_util.c */
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
#define taosolveapplication_ PTAOSOLVEAPPLICATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosolveapplication_ ptaosolveapplication
#else
#define taosolveapplication_ ptaosolveapplication_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosolveapplication_ TAOSOLVEAPPLICATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosolveapplication_ taosolveapplication
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetupapplicationsolver_ PTAOSETUPAPPLICATIONSOLVER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetupapplicationsolver_ ptaosetupapplicationsolver
#else
#define taosetupapplicationsolver_ ptaosetupapplicationsolver_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetupapplicationsolver_ TAOSETUPAPPLICATIONSOLVER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetupapplicationsolver_ taosetupapplicationsolver
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetoptions_ PTAOSETOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetoptions_ ptaosetoptions
#else
#define taosetoptions_ ptaosetoptions_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetoptions_ TAOSETOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetoptions_ taosetoptions
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taocopydualsofvariablebounds_ PTAOCOPYDUALSOFVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taocopydualsofvariablebounds_ ptaocopydualsofvariablebounds
#else
#define taocopydualsofvariablebounds_ ptaocopydualsofvariablebounds_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taocopydualsofvariablebounds_ TAOCOPYDUALSOFVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taocopydualsofvariablebounds_ taocopydualsofvariablebounds
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taosolveapplication_(TAO_APPLICATION *taoapp,TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoSolveApplication(*taoapp,*tao);
}
void PETSC_STDCALL taosetupapplicationsolver_(TAO_APPLICATION *myapp,TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoSetupApplicationSolver(*myapp,*tao);
}
void PETSC_STDCALL taosetoptions_(TAO_APPLICATION *taoapp,TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoSetOptions(*taoapp,*tao);
}
void PETSC_STDCALL taocopydualsofvariablebounds_(TAO_SOLVER *tao,Vec DXL,Vec DXU, int *__ierr ){
*__ierr = TaoCopyDualsOfVariableBounds(*tao,
	(Vec)PetscToPointer((DXL) ),
	(Vec)PetscToPointer((DXU) ));
}
#if defined(__cplusplus)
}
#endif
