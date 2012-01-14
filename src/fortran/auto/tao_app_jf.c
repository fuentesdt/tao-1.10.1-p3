#include "petsc.h"
#include "petscfix.h"
/* tao_app_j.c */
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
#define taoappsetfunctionvec_ PTAOAPPSETFUNCTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetfunctionvec_ ptaoappsetfunctionvec
#else
#define taoappsetfunctionvec_ ptaoappsetfunctionvec_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetfunctionvec_ TAOAPPSETFUNCTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetfunctionvec_ taoappsetfunctionvec
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetjacobianmat_ PTAOAPPSETJACOBIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetjacobianmat_ ptaoappsetjacobianmat
#else
#define taoappsetjacobianmat_ ptaoappsetjacobianmat_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetjacobianmat_ TAOAPPSETJACOBIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetjacobianmat_ taoappsetjacobianmat
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappgetjacobianmat_ PTAOAPPGETJACOBIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappgetjacobianmat_ ptaoappgetjacobianmat
#else
#define taoappgetjacobianmat_ ptaoappgetjacobianmat_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappgetjacobianmat_ TAOAPPGETJACOBIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappgetjacobianmat_ taoappgetjacobianmat
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoappsetfunctionvec_(TAO_APPLICATION *taoapp,Vec r, int *__ierr ){
*__ierr = TaoAppSetFunctionVec(*taoapp,
	(Vec)PetscToPointer((r) ));
}
void PETSC_STDCALL taoappsetjacobianmat_(TAO_APPLICATION *taoapp,Mat J,Mat JP, int *__ierr ){
*__ierr = TaoAppSetJacobianMat(*taoapp,
	(Mat)PetscToPointer((J) ),
	(Mat)PetscToPointer((JP) ));
}
void PETSC_STDCALL taoappgetjacobianmat_(TAO_APPLICATION *taoapp,Mat *J,Mat *JP, int *__ierr ){
*__ierr = TaoAppGetJacobianMat(*taoapp,J,JP);
}
#if defined(__cplusplus)
}
#endif
