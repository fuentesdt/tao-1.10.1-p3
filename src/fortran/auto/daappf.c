#include "petsc.h"
#include "petscfix.h"
/* daapp.c */
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

#include "taodaapplication.h"
#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetdaapp_ PTAOAPPSETDAAPP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetdaapp_ ptaoappsetdaapp
#else
#define taoappsetdaapp_ ptaoappsetdaapp_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetdaapp_ TAOAPPSETDAAPP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetdaapp_ taoappsetdaapp
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetda_ PDAAPPGETDA
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetda_ pdaappgetda
#else
#define daappgetda_ pdaappgetda_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetda_ DAAPPGETDA
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetda_ daappgetda
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetnumberofdagrids_ PDAAPPGETNUMBEROFDAGRIDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetnumberofdagrids_ pdaappgetnumberofdagrids
#else
#define daappgetnumberofdagrids_ pdaappgetnumberofdagrids_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetnumberofdagrids_ DAAPPGETNUMBEROFDAGRIDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetnumberofdagrids_ daappgetnumberofdagrids
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetcurrentlevel_ PDAAPPGETCURRENTLEVEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetcurrentlevel_ pdaappgetcurrentlevel
#else
#define daappgetcurrentlevel_ pdaappgetcurrentlevel_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetcurrentlevel_ DAAPPGETCURRENTLEVEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetcurrentlevel_ daappgetcurrentlevel
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsethessianmat_ PDAAPPSETHESSIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsethessianmat_ pdaappsethessianmat
#else
#define daappsethessianmat_ pdaappsethessianmat_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsethessianmat_ DAAPPSETHESSIANMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsethessianmat_ daappsethessianmat
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetsolution_ PDAAPPGETSOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetsolution_ pdaappgetsolution
#else
#define daappgetsolution_ pdaappgetsolution_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetsolution_ DAAPPGETSOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetsolution_ daappgetsolution
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetinterpolationmatrix_ PDAAPPGETINTERPOLATIONMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetinterpolationmatrix_ pdaappgetinterpolationmatrix
#else
#define daappgetinterpolationmatrix_ pdaappgetinterpolationmatrix_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetinterpolationmatrix_ DAAPPGETINTERPOLATIONMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetinterpolationmatrix_ daappgetinterpolationmatrix
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetvariablebounds_ PDAAPPGETVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetvariablebounds_ pdaappgetvariablebounds
#else
#define daappgetvariablebounds_ pdaappgetvariablebounds_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappgetvariablebounds_ DAAPPGETVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappgetvariablebounds_ daappgetvariablebounds
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetinitialsolution_ PDAAPPSETINITIALSOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetinitialsolution_ pdaappsetinitialsolution
#else
#define daappsetinitialsolution_ pdaappsetinitialsolution_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetinitialsolution_ DAAPPSETINITIALSOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetinitialsolution_ daappsetinitialsolution
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetmattype_ PDAAPPSETMATTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetmattype_ pdaappsetmattype
#else
#define daappsetmattype_ pdaappsetmattype_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetmattype_ DAAPPSETMATTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetmattype_ daappsetmattype
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetoptions_ PDAAPPSETOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetoptions_ pdaappsetoptions
#else
#define daappsetoptions_ pdaappsetoptions_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define daappsetoptions_ DAAPPSETOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define daappsetoptions_ daappsetoptions
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoappsetdaapp_(TAO_APPLICATION *daapplication,DA* tda,PetscInt *nnda, int *__ierr ){
*__ierr = TaoAppSetDAApp(*daapplication,tda,*nnda);
}
void PETSC_STDCALL daappgetda_(TAO_APPLICATION *daapplication,PetscInt *n,DA *da, int *__ierr ){
*__ierr = DAAppGetDA(*daapplication,*n,da);
}
void PETSC_STDCALL daappgetnumberofdagrids_(TAO_APPLICATION *daapplication,PetscInt *n, int *__ierr ){
*__ierr = DAAppGetNumberOfDAGrids(*daapplication,n);
}
void PETSC_STDCALL daappgetcurrentlevel_(TAO_APPLICATION *daapplication,PetscInt *n, int *__ierr ){
*__ierr = DAAppGetCurrentLevel(*daapplication,n);
}
void PETSC_STDCALL daappsethessianmat_(TAO_APPLICATION *daapplication, int *__ierr ){
*__ierr = DAAppSetHessianMat(*daapplication);
}
void PETSC_STDCALL daappgetsolution_(TAO_APPLICATION *daapplication,PetscInt *level,Vec *X, int *__ierr ){
*__ierr = DAAppGetSolution(*daapplication,*level,X);
}
void PETSC_STDCALL daappgetinterpolationmatrix_(TAO_APPLICATION *daapplication,PetscInt *level,Mat *Interpolate,Vec *CScale, int *__ierr ){
*__ierr = DAAppGetInterpolationMatrix(*daapplication,*level,Interpolate,CScale);
}
void PETSC_STDCALL daappgetvariablebounds_(TAO_APPLICATION *daapplication,PetscInt *level,Vec *XL,Vec *XU, int *__ierr ){
*__ierr = DAAppGetVariableBounds(*daapplication,*level,XL,XU);
}
void PETSC_STDCALL daappsetinitialsolution_(TAO_APPLICATION *daapplication,Vec X0, int *__ierr ){
*__ierr = DAAppSetInitialSolution(*daapplication,
	(Vec)PetscToPointer((X0) ));
}
void PETSC_STDCALL daappsetmattype_(TAO_APPLICATION *daapplication, MatType *mattype, int *__ierr ){
*__ierr = DAAppSetMatType(*daapplication,*mattype);
}
void PETSC_STDCALL daappsetoptions_(TAO_APPLICATION *daapplication, int *__ierr ){
*__ierr = DAAppSetOptions(*daapplication);
}
#if defined(__cplusplus)
}
#endif
