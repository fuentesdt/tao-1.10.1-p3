#include "private/fortranimpl.h"
#include "tao_solver.h"

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define taoprintstatement_         TAOPRINTSTATEMENT

#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define taoprintstatement_         taoprintstatement

#endif

EXTERN_C_BEGIN
void PETSC_STDCALL taoprintstatement_(TAO_SOLVER *tao, CHAR statement PETSC_MIXED_LEN(len), int *ierr PETSC_END_LEN(len)) 
{
  char *c1;

  FIXCHAR(statement,len,c1);
  *ierr = TaoPrintStatement(*tao,c1);
  FREECHAR(statement,c1);
}
EXTERN_C_END
