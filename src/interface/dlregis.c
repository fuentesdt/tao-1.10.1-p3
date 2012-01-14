/*$Id$*/

#include "tao_solver.h"

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PetscDLLibraryRegister_tao"
/*
  DLLibraryRegister - This function is called when the dynamic library it is in is opened.

  This registers all of the TAO methods that are in the basic libtao library.

  Input Parameter:
  path - library path
 */
int PetscDLLibraryRegister_tao(const char *path)
{
  int info;
  TaoFunctionBegin;

#ifdef TAO_USE_PETSC
  info = PetscInitializeNoArguments(); if (info) return 1;
#endif

  /*
      If we got here then PETSc was properly loaded
  */
  info = TaoRegisterAll(path);CHKERRQ(info);
  TaoFunctionReturn(0);
}
EXTERN_C_END

