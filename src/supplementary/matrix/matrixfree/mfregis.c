#ifdef PETSC_RCS_HEADER
static char vcid[] = "$Id$";
#endif

#include "src/supplementary/matrix/matrixfree/tao_mfj.h"   /*I  "tao_solver.h"   I*/

EXTERN_C_BEGIN
extern int MatTaoMFCreate_Default(MatTaoMFCtx);
EXTERN_C_END

#undef __FUNC__  
#define __FUNC__ "MatTaoMFRegisterAll"
/*@C
  MatTaoMFRegisterAll - Registers all of the differencing interval computation
  routines in the MatTaoMF package.

  Not Collective

  Level: developer

.keywords: MatTaoMF, register, all

.seealso:  MatTaoMFRegisterDestroy(), MatTaoMFRegister(), MatTaoMFCreate(), 
           MatTaoMFSetType()
@*/
int MatTaoMFRegisterAll(char *path)
{
  int info;

  PetscFunctionBegin;
  MatTaoMFRegisterAllCalled = 1;

  info = MatTaoMFRegisterDynamic("default",path,"MatTaoMFCreate_Default",MatTaoMFCreate_Default);CHKERRQ(info);
  PetscFunctionReturn(0);
}


