#include "petsc.h"
#include "petscfix.h"
/* taoapp_ga.c */
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

#include "taoapp_ga.h"
#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetupgaapplicationsolver_ PTAOSETUPGAAPPLICATIONSOLVER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetupgaapplicationsolver_ ptaosetupgaapplicationsolver
#else
#define taosetupgaapplicationsolver_ ptaosetupgaapplicationsolver_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetupgaapplicationsolver_ TAOSETUPGAAPPLICATIONSOLVER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetupgaapplicationsolver_ taosetupgaapplicationsolver
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosolvegaapplication_ PTAOSOLVEGAAPPLICATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosolvegaapplication_ ptaosolvegaapplication
#else
#define taosolvegaapplication_ ptaosolvegaapplication_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosolvegaapplication_ TAOSOLVEGAAPPLICATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosolvegaapplication_ taosolvegaapplication
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogaappsetinitialsolutionvec_ PTAOGAAPPSETINITIALSOLUTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogaappsetinitialsolutionvec_ ptaogaappsetinitialsolutionvec
#else
#define taogaappsetinitialsolutionvec_ ptaogaappsetinitialsolutionvec_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogaappsetinitialsolutionvec_ TAOGAAPPSETINITIALSOLUTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogaappsetinitialsolutionvec_ taogaappsetinitialsolutionvec
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogaappmonitor_ PTAOGAAPPMONITOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogaappmonitor_ ptaogaappmonitor
#else
#define taogaappmonitor_ ptaogaappmonitor_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taogaappmonitor_ TAOGAAPPMONITOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taogaappmonitor_ taogaappmonitor
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taosetupgaapplicationsolver_(TAO_GA_APPLICATION *myapp,TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoSetupGAApplicationSolver(*myapp,*tao);
}
void PETSC_STDCALL taosolvegaapplication_(TAO_GA_APPLICATION *gaapp,TAO_SOLVER *tao, int *__ierr ){
*__ierr = TaoSolveGAApplication(*gaapp,*tao);
}
void PETSC_STDCALL taogaappsetinitialsolutionvec_(TAO_GA_APPLICATION *gaapp,GAVec *xx, int *__ierr ){
*__ierr = TaoGAAppSetInitialSolutionVec(*gaapp,*xx);
}
void PETSC_STDCALL taogaappmonitor_(TAO_GA_APPLICATION *gaapp, int *__ierr ){
*__ierr = TaoGAAppMonitor(*gaapp);
}
#if defined(__cplusplus)
}
#endif
