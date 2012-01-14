/*$Id$*/

#include "private/fortranimpl.h" 
/* #include "sys.h" */
#include "petscsys.h"

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoinitializefortran_       TAOINITIALIZEFORTRAN
#define taosetcommonblock_          TAOSETCOMMONBLOCK
#define tao_null_function_          TAO_NULL_FUNCTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define taoinitializefortran_       taoinitializefortran
#define taosetcommonblock_          taosetcommonblock
#define tao_null_function_          tao_null_function
#endif

#if defined(PETSC_HAVE_FORTRAN_UNDERSCORE_UNDERSCORE)
#define tao_null_function_  tao_null_function__
#endif

EXTERN_C_BEGIN
extern void PETSC_STDCALL taosetcommonblock_(void);
EXTERN_C_END

/*@C
   TaoInitializeFortran - Routine that should be called from C after
   the call to TaoInitialize() if one is using a C main program
   that calls Fortran routines that in turn call TAO routines.

   Collective on MPI_COMM_WORLD

   Level: intermediate

   Notes:
   TaoInitializeFortran() initializes some of the default TAO variables
   for use in Fortran if a user's main program is written in C.  
   TaoInitializeFortran() is NOT needed if a user's main
   program is written in Fortran; in this case, just calling
   TaoInitialize() in the main (Fortran) program is sufficient.

.seealso:  TaoInitialize()

.keywords: Mixing C and Fortran, passing TAO objects to Fortran
@*/

int TaoInitializeFortran(void)
{
  taosetcommonblock_();
  return 0;
}
  
EXTERN_C_BEGIN

void PETSC_STDCALL taoinitializefortran_(int *info)
{
  *info = TaoInitializeFortran();
}

/*
  A valid address for the Fortran variable TAO_NULL_FUNCTION
*/
void tao_null_function_(void)
{
  return;
}

EXTERN_C_END

