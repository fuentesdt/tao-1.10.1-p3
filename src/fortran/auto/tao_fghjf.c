#include "petsc.h"
#include "petscfix.h"
/* tao_fghj.c */
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
#define taoincrementgradientscounter_ PTAOINCREMENTGRADIENTSCOUNTER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoincrementgradientscounter_ ptaoincrementgradientscounter
#else
#define taoincrementgradientscounter_ ptaoincrementgradientscounter_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoincrementgradientscounter_ TAOINCREMENTGRADIENTSCOUNTER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoincrementgradientscounter_ taoincrementgradientscounter
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoevaluatevariablebounds_ PTAOEVALUATEVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoevaluatevariablebounds_ ptaoevaluatevariablebounds
#else
#define taoevaluatevariablebounds_ ptaoevaluatevariablebounds_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoevaluatevariablebounds_ TAOEVALUATEVARIABLEBOUNDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoevaluatevariablebounds_ taoevaluatevariablebounds
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetlagrangiangradientvector_ PTAOSETLAGRANGIANGRADIENTVECTOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetlagrangiangradientvector_ ptaosetlagrangiangradientvector
#else
#define taosetlagrangiangradientvector_ ptaosetlagrangiangradientvector_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taosetlagrangiangradientvector_ TAOSETLAGRANGIANGRADIENTVECTOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taosetlagrangiangradientvector_ taosetlagrangiangradientvector
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoincrementgradientscounter_(TAO_SOLVER *solver,TaoInt *nevals, int *__ierr ){
*__ierr = TaoIncrementGradientsCounter(*solver,*nevals);
}
void PETSC_STDCALL taoevaluatevariablebounds_(TAO_SOLVER *tao,TaoVec *xxll,TaoVec *xxuu, int *__ierr ){
*__ierr = TaoEvaluateVariableBounds(*tao,
	(TaoVec* )PetscToPointer((xxll) ),
	(TaoVec* )PetscToPointer((xxuu) ));
}
void PETSC_STDCALL taosetlagrangiangradientvector_(TAO_SOLVER *solver,TaoVec* gg, int *__ierr ){
*__ierr = TaoSetLagrangianGradientVector(*solver,
	(TaoVec* )PetscToPointer((gg) ));
}
#if defined(__cplusplus)
}
#endif
