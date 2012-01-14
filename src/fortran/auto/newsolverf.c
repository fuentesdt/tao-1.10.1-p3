#include "petsc.h"
#include "petscfix.h"
/* newsolver.c */
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
#define taosetmethodfromoptions_ PTAOSETMETHODFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetmethodfromoptions_ ptaosetmethodfromoptions
#else
#define taosetmethodfromoptions_ ptaosetmethodfromoptions_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetmethodfromoptions_ TAOSETMETHODFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetmethodfromoptions_ taosetmethodfromoptions
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetfromoptions_ PTAOSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetfromoptions_ ptaosetfromoptions
#else
#define taosetfromoptions_ ptaosetfromoptions_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetfromoptions_ TAOSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetfromoptions_ taosetfromoptions
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdefaultparameters_ PTAOSETDEFAULTPARAMETERS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdefaultparameters_ ptaosetdefaultparameters
#else
#define taosetdefaultparameters_ ptaosetdefaultparameters_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdefaultparameters_ TAOSETDEFAULTPARAMETERS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdefaultparameters_ taosetdefaultparameters
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdefaultstatistics_ PTAOSETDEFAULTSTATISTICS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdefaultstatistics_ ptaosetdefaultstatistics
#else
#define taosetdefaultstatistics_ ptaosetdefaultstatistics_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdefaultstatistics_ TAOSETDEFAULTSTATISTICS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdefaultstatistics_ taosetdefaultstatistics
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdefaultmonitors_ PTAOSETDEFAULTMONITORS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdefaultmonitors_ ptaosetdefaultmonitors
#else
#define taosetdefaultmonitors_ ptaosetdefaultmonitors_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdefaultmonitors_ TAOSETDEFAULTMONITORS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdefaultmonitors_ taosetdefaultmonitors
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetup_ PTAOSETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetup_ ptaosetup
#else
#define taosetup_ ptaosetup_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetup_ TAOSETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetup_ taosetup
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taodestroy_ PTAODESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taodestroy_ ptaodestroy
#else
#define taodestroy_ ptaodestroy_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taodestroy_ TAODESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taodestroy_ taodestroy
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdown_ PTAOSETDOWN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdown_ ptaosetdown
#else
#define taosetdown_ ptaosetdown_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetdown_ TAOSETDOWN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetdown_ taosetdown
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoresetsolver_ PTAORESETSOLVER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoresetsolver_ ptaoresetsolver
#else
#define taoresetsolver_ ptaoresetsolver_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoresetsolver_ TAORESETSOLVER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoresetsolver_ taoresetsolver
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taosetmethodfromoptions_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoSetMethodFromOptions(*solver);
}
void PETSC_STDCALL taosetfromoptions_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoSetFromOptions(*solver);
}
void PETSC_STDCALL taosetdefaultparameters_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoSetDefaultParameters(*solver);
}
void PETSC_STDCALL taosetdefaultstatistics_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoSetDefaultStatistics(*solver);
}
void PETSC_STDCALL taosetdefaultmonitors_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoSetDefaultMonitors(*solver);
}
void PETSC_STDCALL taosetup_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoSetUp(*solver);
}
void PETSC_STDCALL taodestroy_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoDestroy(*solver);
}
void PETSC_STDCALL taosetdown_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoSetDown(*solver);
}
void PETSC_STDCALL taoresetsolver_(TAO_SOLVER *solver, int *__ierr ){
*__ierr = TaoResetSolver(*solver);
}
#if defined(__cplusplus)
}
#endif
