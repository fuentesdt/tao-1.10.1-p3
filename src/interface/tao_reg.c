/*$Id$*/

#include "src/tao_impl.h"     /*I  "tao_solver.h"  I*/

EXTERN_C_BEGIN
extern int TaoCreate_LMVM(TAO_SOLVER);
extern int TaoCreate_NLS(TAO_SOLVER);
extern int TaoCreate_NTR(TAO_SOLVER);
extern int TaoCreate_NTL(TAO_SOLVER);
extern int TaoCreate_CG(TAO_SOLVER);

extern int TaoCreate_TRON(TAO_SOLVER);
extern int TaoCreate_BQPIP(TAO_SOLVER);
extern int TaoCreate_BLMVM(TAO_SOLVER);
extern int TaoCreate_BNLS(TAO_SOLVER);
extern int TaoCreate_GPCG(TAO_SOLVER);
extern int TaoCreate_QPIP(TAO_SOLVER);

extern int TaoCreate_NLSQ(TAO_SOLVER);
extern int TaoCreate_BLM(TAO_SOLVER);
extern int TaoCreate_SSILS(TAO_SOLVER);
extern int TaoCreate_SSFLS(TAO_SOLVER);
extern int TaoCreate_ASILS(TAO_SOLVER);
extern int TaoCreate_ASFLS(TAO_SOLVER);
extern int TaoCreate_ISILS(TAO_SOLVER);
extern int TaoCreate_KT(TAO_SOLVER);
extern int TaoCreate_BCG(TAO_SOLVER);
extern int TaoCreate_RSCS(TAO_SOLVER);
extern int TaoCreate_ICP(TAO_SOLVER);
extern int TaoCreate_NelderMead(TAO_SOLVER);
extern int TaoCreate_FD(TAO_SOLVER);
EXTERN_C_END

/* #undef USE_DYNAMIC_LIBRARIES */


/*
    This routine is used by TaoSetType() to make sure that 
    TaoRegisterAll() is called at least once. In general, if 
    there is more than one DLL, then TaoRegisterAll() may be
    called several times.
*/
static int TaoRegisterAllCalled=0;

int TaoStandardRegisterAll(){
  int info;
  char path[PETSC_MAX_PATH_LEN];
  info = PetscStrcpy(path,TAO_LIB_DIR); CHKERRQ(info);
  info = PetscStrcat(path,"/libtao");
  TaoFunctionBegin;
  if (TaoRegisterAllCalled){
    TaoFunctionReturn(0);
  }
#ifdef PETSC_USE_DYNAMIC_LIBRARIES
  info = PetscDLLibraryAppend(PETSC_COMM_WORLD,&DLLibrariesLoaded,path);CHKERRQ(info);
#endif
  info=TaoRegisterAll(path); CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoRegisterAll"
/*@C
   TaoRegisterAll - Registers all of the minimization methods in the TAO_SOLVER package.

   Not Collective

   Level: advanced

.keywords: TAO_SOLVER, register, all

.seealso:  TaoRegisterDestroy()
@*/
int TaoRegisterAll(const char *path)
{
  int info;
  TaoFunctionBegin;
  TaoRegisterAllCalled = 1;

  info = TaoRegisterDynamic("tao_lmvm",path,"TaoCreate_LMVM",TaoCreate_LMVM); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_nls",path,"TaoCreate_NLS",TaoCreate_NLS); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_cg",path,"TaoCreate_CG",TaoCreate_CG); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_bqpip",path,"TaoCreate_BQPIP",TaoCreate_BQPIP); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_blmvm",path,"TaoCreate_BLMVM",TaoCreate_BLMVM); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_bnls",path,"TaoCreate_BNLS",TaoCreate_BNLS); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_tron",path,"TaoCreate_TRON",TaoCreate_TRON); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_nm",path,"TaoCreate_NelderMead",TaoCreate_NelderMead); CHKERRQ(info);
#ifdef TAO_USE_PETSC
  info = TaoRegisterDynamic("tao_ntl",path,"TaoCreate_NTL",TaoCreate_NTL); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_ntr",path,"TaoCreate_NTR",TaoCreate_NTR); CHKERRQ(info);
#endif
  info = TaoRegisterDynamic("tao_gpcg",path,"TaoCreate_GPCG",TaoCreate_GPCG); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_nm",path,"TaoCreate_NelderMead",TaoCreate_NelderMead); CHKERRQ(info);

  //  info = TaoRegisterDynamic("tao_nlsq",path,"TaoCreate_NLSQ",TaoCreate_NLSQ); CHKERRQ(info);

  /* Add registration for the semismooth code using a linesearch. */
  info = TaoRegisterDynamic("tao_ssils",path,"TaoCreate_SSILS",TaoCreate_SSILS); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_ssfls",path,"TaoCreate_SSFLS",TaoCreate_SSFLS); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_asils",path,"TaoCreate_ASILS",TaoCreate_ASILS); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_asfls",path,"TaoCreate_ASFLS",TaoCreate_ASFLS); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_isils",path,"TaoCreate_ISILS",TaoCreate_ISILS); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_kt",path,"TaoCreate_KT",TaoCreate_KT); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_rscs",path,"TaoCreate_RSCS",TaoCreate_RSCS); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_icp",path,"TaoCreate_ICP",TaoCreate_ICP); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_fd_test",path,"TaoCreate_FD",TaoCreate_FD); CHKERRQ(info);
  /*
  info = TaoRegisterDynamic("tao_bcg",path,"TaoCreate_BCG",TaoCreate_BCG); CHKERRQ(info);
  info = TaoRegisterDynamic("tao_qpip",path,"TaoCreate_QPIP",TaoCreate_QPIP); CHKERRQ(info);
  */
  
  TaoFunctionReturn(0);
}

