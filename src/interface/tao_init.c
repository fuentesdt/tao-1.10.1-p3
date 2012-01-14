/*$Id$*/

#include "tao_solver.h"   /*I  "tao_solver.h" I*/

int TaoRegisterEvents();

/* ------------------------Nasty global variables -------------------------------*/
int TaoInitializeCalled = 0;
int TAO_COOKIE = 0;

static int  TaoGlobalArgc   = 0;
static char **TaoGlobalArgs = 0;

#undef __FUNCT__  
#define __FUNCT__ "TaoInitialize"
/*@C 
  TaoInitialize - Initializes the TAO component and many of the packages associated with it.

   Collective on MPI_COMM_WORLD

   Input Parameters:
+  argc - [optional] count of number of command line arguments
.  args - [optional] the command line arguments
.  file - [optional] PETSc database file, defaults to ~username/.petscrc
          (use TAO_NULL for default)
-  help - [optional] Help message to print, use TAO_NULL for no message

   Note:
   TaoInitialize() should always be called near the beginning of your 
   program.  However, this command should come after PetscInitialize()

   Note:
   The input arguments are required if the options database is to be used.

   Level: beginner

.keywords: TAO_SOLVER, initialize

.seealso: TaoInitializeFortran(), TaoFinalize(), PetscInitialize()
@*/
int TaoInitialize(int *argc,char ***args,char file[],const char help[])
{
  int info=0;

  TaoFunctionBegin;

  if (TaoInitializeCalled){ TaoFunctionReturn(0);}
  TaoInitializeCalled++;

  if (argc && args){
    TaoGlobalArgc = *argc;
    TaoGlobalArgs = *args;
  }

  TAO_COOKIE = 0;
  info=TaoLogClassRegister(&TAO_COOKIE,"TAO Solver"); CHKERRQ(info);

  info = TaoRegisterEvents(); CHKERRQ(info);
  info = TaoStandardRegisterAll();CHKERRQ(info);
  info = PetscInfo(0,"TaoInitialize:TAO successfully started\n"); CHKERRQ(info);
  TaoFunctionReturn(info);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoFinalize"
/*@
   TaoFinalize - Checks for options at the end of the TAO program
   and finalizes interfaces with other packages.

   Collective on MPI_COMM_WORLD

   Level: beginner

.keywords: finalize, exit, end

.seealso: TaoInitialize(), PetscFinalize()
@*/
int TaoFinalize(void)
{
  int info;
  
  TaoFunctionBegin;
  TaoInitializeCalled--;
  if (TaoInitializeCalled==0){
    info = PetscInfo(0,"TaoFinalize:Tao successfully ended!\n"); 
           CHKERRQ(info);
    info = TaoRegisterDestroy(); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoRegisterEvents"
// int Tao_Solve, Tao_LineSearch;
/*
   TaoRegisterEvents - Registers TAO events for use in performance logging.
*/
int TaoRegisterEvents()
{
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoGetArgs"
/*@C
   TaoGetArgs - Allows you to access the raw command line arguments anywhere
     after TaoInitialize() is called but before TaoFinalize().

   Not Collective

   Output Parameters:
+  argc - count of number of command line arguments
-  args - the command line arguments

   Level: developer

   Notes:
      This is usually used to pass the command line arguments into other 
      libraries that are called internally deep in TAO or the application.

.keywords: command line arguments
   
.seealso: TaoFinalize(), TaoInitializeFortran()

@*/
int TaoGetArgs(int *argc,char ***args)
{
  TaoFunctionBegin;
  if (!TaoGlobalArgs) {
    SETERRQ(1,"You must call after TaoInitialize()");
  }
  if (argc && args){
    *argc = TaoGlobalArgc;
    *args = TaoGlobalArgs;
  }
  TaoFunctionReturn(0);
}


