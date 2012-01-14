#include "petsc.h"
#include "petscfix.h"
/* tao_mfj.c */
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
#define mattaomfsetfromoptions_ PMATTAOMFSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfsetfromoptions_ pmattaomfsetfromoptions
#else
#define mattaomfsetfromoptions_ pmattaomfsetfromoptions_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfsetfromoptions_ MATTAOMFSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfsetfromoptions_ mattaomfsetfromoptions
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfgeth_ PMATTAOMFGETH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfgeth_ pmattaomfgeth
#else
#define mattaomfgeth_ pmattaomfgeth_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfgeth_ MATTAOMFGETH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfgeth_ mattaomfgeth
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfsetgradienterror_ PMATTAOMFSETGRADIENTERROR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfsetgradienterror_ pmattaomfsetgradienterror
#else
#define mattaomfsetgradienterror_ pmattaomfsetgradienterror_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfsetgradienterror_ MATTAOMFSETGRADIENTERROR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfsetgradienterror_ mattaomfsetgradienterror
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfsethhistory_ PMATTAOMFSETHHISTORY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfsethhistory_ pmattaomfsethhistory
#else
#define mattaomfsethhistory_ pmattaomfsethhistory_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfsethhistory_ MATTAOMFSETHHISTORY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfsethhistory_ mattaomfsethhistory
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfresethhistory_ PMATTAOMFRESETHHISTORY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfresethhistory_ pmattaomfresethhistory
#else
#define mattaomfresethhistory_ pmattaomfresethhistory_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mattaomfresethhistory_ MATTAOMFRESETHHISTORY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mattaomfresethhistory_ mattaomfresethhistory
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL mattaomfsetfromoptions_(Mat mat, int *__ierr ){
*__ierr = MatTaoMFSetFromOptions(
	(Mat)PetscToPointer((mat) ));
}
void PETSC_STDCALL mattaomfgeth_(Mat mat,Scalar *h, int *__ierr ){
*__ierr = MatTaoMFGetH(
	(Mat)PetscToPointer((mat) ),
	(Scalar* )PetscToPointer((h) ));
}
void PETSC_STDCALL mattaomfsetgradienterror_(Mat mat,double *error, int *__ierr ){
*__ierr = MatTaoMFSetGradientError(
	(Mat)PetscToPointer((mat) ),*error);
}
void PETSC_STDCALL mattaomfsethhistory_(Mat J,Scalar *history,int *nhistory, int *__ierr ){
*__ierr = MatTaoMFSetHHistory(
	(Mat)PetscToPointer((J) ),
	(Scalar* )PetscToPointer((history) ),*nhistory);
}
void PETSC_STDCALL mattaomfresethhistory_(Mat J, int *__ierr ){
*__ierr = MatTaoMFResetHHistory(
	(Mat)PetscToPointer((J) ));
}
#if defined(__cplusplus)
}
#endif
