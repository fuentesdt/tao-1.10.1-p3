#include "petsc.h"
#include "petscfix.h"
/* tao_app.c */
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
#define taoappdestroy_ PTAOAPPDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappdestroy_ ptaoappdestroy
#else
#define taoappdestroy_ ptaoappdestroy_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappdestroy_ TAOAPPDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappdestroy_ taoappdestroy
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetinitialsolutionvec_ PTAOAPPSETINITIALSOLUTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetinitialsolutionvec_ ptaoappsetinitialsolutionvec
#else
#define taoappsetinitialsolutionvec_ ptaoappsetinitialsolutionvec_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetinitialsolutionvec_ TAOAPPSETINITIALSOLUTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetinitialsolutionvec_ taoappsetinitialsolutionvec
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetdefaultsolutionvec_ PTAOAPPSETDEFAULTSOLUTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetdefaultsolutionvec_ ptaoappsetdefaultsolutionvec
#else
#define taoappsetdefaultsolutionvec_ ptaoappsetdefaultsolutionvec_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetdefaultsolutionvec_ TAOAPPSETDEFAULTSOLUTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetdefaultsolutionvec_ taoappsetdefaultsolutionvec
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappgetsolutionvec_ PTAOAPPGETSOLUTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappgetsolutionvec_ ptaoappgetsolutionvec
#else
#define taoappgetsolutionvec_ ptaoappgetsolutionvec_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappgetsolutionvec_ TAOAPPGETSOLUTIONVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappgetsolutionvec_ taoappgetsolutionvec
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappmonitor_ PTAOAPPMONITOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappmonitor_ ptaoappmonitor
#else
#define taoappmonitor_ ptaoappmonitor_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappmonitor_ TAOAPPMONITOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappmonitor_ taoappmonitor
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappaddobject_ PTAOAPPADDOBJECT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappaddobject_ ptaoappaddobject
#else
#define taoappaddobject_ ptaoappaddobject_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappaddobject_ TAOAPPADDOBJECT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappaddobject_ taoappaddobject
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappqueryforobject_ PTAOAPPQUERYFOROBJECT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappqueryforobject_ ptaoappqueryforobject
#else
#define taoappqueryforobject_ ptaoappqueryforobject_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappqueryforobject_ TAOAPPQUERYFOROBJECT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappqueryforobject_ taoappqueryforobject
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappqueryremoveobject_ PTAOAPPQUERYREMOVEOBJECT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappqueryremoveobject_ ptaoappqueryremoveobject
#else
#define taoappqueryremoveobject_ ptaoappqueryremoveobject_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappqueryremoveobject_ TAOAPPQUERYREMOVEOBJECT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappqueryremoveobject_ taoappqueryremoveobject
#endif
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetfromoptions_ PTAOAPPSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetfromoptions_ ptaoappsetfromoptions
#else
#define taoappsetfromoptions_ ptaoappsetfromoptions_
#endif
#else
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoappsetfromoptions_ TAOAPPSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define taoappsetfromoptions_ taoappsetfromoptions
#endif
#endif



/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL taoappdestroy_(TAO_APPLICATION *taoapp, int *__ierr ){
*__ierr = TaoAppDestroy(*taoapp);
}
void PETSC_STDCALL taoappsetinitialsolutionvec_(TAO_APPLICATION *taoapp,Vec xx, int *__ierr ){
*__ierr = TaoAppSetInitialSolutionVec(*taoapp,
	(Vec)PetscToPointer((xx) ));
}
void PETSC_STDCALL taoappsetdefaultsolutionvec_(TAO_APPLICATION *taoapp,Vec xx, int *__ierr ){
*__ierr = TaoAppSetDefaultSolutionVec(*taoapp,
	(Vec)PetscToPointer((xx) ));
}
void PETSC_STDCALL taoappgetsolutionvec_(TAO_APPLICATION *taoapp,Vec *X, int *__ierr ){
*__ierr = TaoAppGetSolutionVec(*taoapp,X);
}
void PETSC_STDCALL taoappmonitor_(TAO_APPLICATION *taoapp, int *__ierr ){
*__ierr = TaoAppMonitor(*taoapp);
}
void PETSC_STDCALL taoappaddobject_(TAO_APPLICATION *taoapp,char *key,void*ctx,TaoInt *id, int *__ierr ){
*__ierr = TaoAppAddObject(*taoapp,key,ctx,
	(TaoInt* )PetscToPointer((id) ));
}
void PETSC_STDCALL taoappqueryforobject_(TAO_APPLICATION *taoapp,char *key,void**ctx, int *__ierr ){
*__ierr = TaoAppQueryForObject(*taoapp,key,ctx);
}
void PETSC_STDCALL taoappqueryremoveobject_(TAO_APPLICATION *taoapp,char *key, int *__ierr ){
*__ierr = TaoAppQueryRemoveObject(*taoapp,key);
}
void PETSC_STDCALL taoappsetfromoptions_(TAO_APPLICATION *taoapp, int *__ierr ){
*__ierr = TaoAppSetFromOptions(*taoapp);
}
#if defined(__cplusplus)
}
#endif
