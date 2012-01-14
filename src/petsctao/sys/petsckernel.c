#include "tao_general.h"  /*I "tao_general.h"  I*/
#ifdef TAO_USE_PETSC
#include "src/tao_impl.h"

TaoTruth TaoBeganPetsc = TAO_FALSE; 

// extern int PetscInitializeCalled=0;
#include "petscsys.h"


static PetscFList TaoList = 0;

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintVersion"
/*
   TaoPrintVersion - Prints TAO version info.

   Collective on MPI_Comm
*/
int TaoPrintVersion(MPI_Comm comm)
{
  int info=0;
  
  PetscFunctionBegin;

  info = PetscHelpPrintf(comm,"--------------------------------------------\
------------------------------\n"); CHKERRQ(info);
  info = PetscHelpPrintf(comm,"\t   %s\n",TAO_VERSION_NUMBER); CHKERRQ(info);
  info = PetscHelpPrintf(comm,"%s",TAO_AUTHOR_INFO); CHKERRQ(info);
  info = PetscHelpPrintf(comm,"See docs/manualpages/index.html for help. \n"); CHKERRQ(info);
#if !defined(PARCH_win32)
  info = PetscHelpPrintf(comm,"TAO libraries linked from %s\n",TAO_LIB_DIR); CHKERRQ(info);
#endif
  info = PetscHelpPrintf(comm,"--------------------------------------------\
------------------------------\n"); CHKERRQ(info);

  PetscFunctionReturn(info);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintHelpIntro"
/*
   TaoPrintHelpIntro - Prints introductory TAO help info.

   Collective on MPI_Comm
*/
int TaoPrintHelpIntro(MPI_Comm comm)
{
  int info=0;
  
  PetscFunctionBegin;

  info = PetscHelpPrintf(comm,"--------------------------------------------\
------------------------------\n"); CHKERRQ(info);
  info = PetscHelpPrintf(comm,"TAO help information includes that for the PETSc libraries, which provide\n"); CHKERRQ(info);
  info = PetscHelpPrintf(comm,"low-level system infrastructure and linear algebra tools.\n"); CHKERRQ(info);
  info = PetscHelpPrintf(comm,"--------------------------------------------\
------------------------------\n"); CHKERRQ(info);

  PetscFunctionReturn(info);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintStatement"
/*@
  TaoPrintStatement - prints a character string to stdout.

  Not Collective
  
  Input Parameters:
+ tao - the TAO_SOLVER solver context
- statement - the string to print

  Level: beginner
@*/
int TaoPrintStatement(TAO_SOLVER tao, const char *statement){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  PetscPrintf(((PetscObject)tao)->comm,statement);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintInt"
int TaoPrintInt(TAO_SOLVER tao, const char *statement, TaoInt n){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  PetscPrintf(((PetscObject)tao)->comm,statement,(PetscInt)n);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintDouble"
int TaoPrintDouble(TAO_SOLVER tao, const char *statement,double dd){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  PetscPrintf(((PetscObject)tao)->comm,statement,dd);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintString"
int TaoPrintString(TAO_SOLVER tao, const char *statement,const char *str){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  PetscPrintf(((PetscObject)tao)->comm,statement,str);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionsHead"
int TaoOptionsHead(const char *heading){
  int info;
  PetscFunctionBegin;
  info = PetscOptionsHead(heading);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionsTail"
int TaoOptionsTail(){
  int info;
  PetscFunctionBegin;
  info = PetscOptionsTail();CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintOptionInt"
int TaoOptionInt(const char *opt,const char *text,const char *man,TaoInt defaultv,TaoInt *value,TaoTruth *set){
  int info;
  PetscTruth flg=PETSC_FALSE;
  PetscFunctionBegin;
  info = PetscOptionsInt(opt,text,man,defaultv,value,&flg);CHKERRQ(info);
  if (set){
    if (flg==PETSC_TRUE) *set=TAO_TRUE; else *set=TAO_FALSE;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintOptionDouble"
int TaoOptionDouble(const char *opt,const char *text,const char *man,double defaultv,double *value,TaoTruth *set){
  int info;
  PetscReal pdv=(PetscReal)defaultv, pv;
  PetscTruth flg=PETSC_FALSE;
  PetscFunctionBegin;
  info = PetscOptionsReal(opt,text,man,pdv,&pv,&flg);CHKERRQ(info);
  if (set){
    if (flg==PETSC_TRUE) *set=TAO_TRUE; else *set=TAO_FALSE;
  }
  if (value&&flg==PETSC_TRUE){
    *value=(double)pv;
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintOptionString"
int TaoOptionString(const char *opt,const char *text,const char *man,const char* defaultv,char *value, TaoInt len, TaoTruth *set){
  int info;
  PetscTruth flg=PETSC_FALSE;
  PetscFunctionBegin;
  info = PetscOptionsString(opt,text,man,defaultv,value,len,&flg);CHKERRQ(info);
  if (set){
    if (flg==PETSC_TRUE) *set=TAO_TRUE; else *set=TAO_FALSE;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintOptionName"
int TaoOptionName(const char *opt,const char *text,const char *man,TaoTruth *set){
  int info;
  PetscTruth flg=PETSC_FALSE;
  PetscFunctionBegin;
  info = PetscOptionsName(opt,text,man,&flg);CHKERRQ(info);
  if (set){
    if (flg==PETSC_TRUE) *set=TAO_TRUE; else *set=TAO_FALSE;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionList"
int TaoOptionList(const char *opt, const char *ltext, const char *man, 
                  const char **list, TaoInt nlist, const char *defaultv, 
                  TaoInt *value, TaoTruth *set)
{
  int info;
  PetscTruth flg=PETSC_FALSE;
  PetscFunctionBegin;
  info = PetscOptionsEList(opt, ltext, man, list, nlist, defaultv, value, &flg); CHKERRQ(info);
  
  if (set) {
    if (PETSC_TRUE == flg) {
      *set=TAO_TRUE; 
    }
    else {
      *set=TAO_FALSE;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoMethodsList"
int TaoMethodsList(const char *opt,const char *ltext,const char *man,const char *defaultv,char *value,TaoInt len,TaoTruth *set){
  int info;
  PetscTruth flg=PETSC_FALSE;
  PetscFunctionBegin;
  info = PetscOptionsList(opt,ltext,man,TaoList,defaultv,value,len,&flg); CHKERRQ(info);
  if (set){
    if (flg==PETSC_TRUE) *set=TAO_TRUE; else *set=TAO_FALSE;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoFindSolver"
int TaoFindSolver(TAO_SOLVER tao, TaoMethod type,  int (**r)(TAO_SOLVER) ){
  int info;
  PetscFunctionBegin;
  info = PetscFListFind(TaoList,((PetscObject)tao)->comm,type,(void (**)(void))r);CHKERRQ(info);
  if (*r){
    info = PetscObjectChangeTypeName((PetscObject)tao,type);CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoRegisterDestroy"
/*@C
   TaoRegisterDestroy - Frees the list of minimization solvers that were
   registered by TaoRegisterDynamic().

   Not Collective

   Level: advanced

.keywords: TAO_SOLVER, register, destroy

.seealso: TaoRegisterAll()
@*/
int TaoRegisterDestroy(){
  int info;
  PetscFunctionBegin;
  if (TaoList) {
    info=PetscFListDestroy(&TaoList);CHKERRQ(info);
    TaoList=0;
  }
  if (TaoBeganPetsc) {
    info = PetscFinalize();CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoRegisterDynamic"
/* --------------------------------------------------------------------- */
/*MC
   TaoRegisterDynamic - Adds a method to the TAO_SOLVER package for unconstrained minimization.

   Synopsis:
   TaoRegisterDynamic(char *name_solver,char *path,char *name_Create,int (*routine_Create)(TAO_SOLVER))

   Not collective

   Input Parameters:
+  name_solver - name of a new user-defined solver
.  path - path (either absolute or relative) the library containing this solver
.  name_Create - name of routine to Create method context
-  routine_Create - routine to Create method context

   Notes:
   TaoRegisterDynamic() may be called multiple times to add several user-defined solvers.

   If dynamic libraries are used, then the fourth input argument (routine_Create)
   is ignored.

   Environmental variables such as ${TAO_DIR}, ${PETSC_ARCH}, ${PETSC_DIR}, 
   and others of the form ${any_environmental_variable} occuring in pathname will be 
   replaced with the appropriate values used when PETSc and TAO were compiled.

   Sample usage:
.vb
   TaoRegisterDynamic("my_solver","/home/username/my_lib/lib/libg_c++/solaris/mylib.a",
                "MySolverCreate",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     TaoSetMethod(solver,"my_solver")
   or at runtime via the option
$     -tao_method my_solver

   Level: advanced

.keywords: TAO_SOLVER, register

.seealso: TaoRegisterAll(), TaoRegisterDestroy()
M*/
#undef __FUNCT__  
#define __FUNCT__ "TaoRegisterDynamic"
EXTERN_C_BEGIN   // make extern c to allow extern c function pointer as argument
int TaoRegisterDynamic(const char *sname,const char *path,const char *name,int (*function)(TAO_SOLVER))
{
  char fullname[256];
  int  info;

  PetscFunctionBegin;
  info = PetscFListConcat(path,name,fullname); CHKERRQ(info);
  //#if defined(PETSC_USE_DYNAMIC_LIBRARIES)
  //info = PetscFListAddDynamic(&TaoList,sname,fullname, 0);CHKERRQ(info);
  //#else
  info = PetscFListAddDynamic(&TaoList,sname,fullname, (void (*)())function);CHKERRQ(info);
  //#endif
  PetscFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "TaoCompareMethod"
/*@C
   TaoCompareMethod - Determines whether the TAO_SOLVER structure is of a
   specified type.

   Not Collective

   Input Parameter:
.  tao - the TAO_SOLVER solver context
.  type - a TAO_SOLVER solver method

   Output Parameter:
.  same - TAO_TRUE if 'tao' is of method 'type'

   Level: developer

.keywords: method
@*/
int TaoCompareMethod(TAO_SOLVER tao, TaoMethod type,TaoTruth *issame){
  int info;
  PetscTruth flag;
  TaoMethod currenttype;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = TaoGetMethod(tao,&currenttype);CHKERRQ(info);
  info = PetscStrcmp(type,currenttype,&flag);CHKERRQ(info);
  if (issame){
    if (flag==PETSC_TRUE) *issame=TAO_TRUE; else *issame=TAO_FALSE;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoStrcmp"
int TaoStrcmp(const char *str1,const char *str2,TaoTruth *flag){
  int info;
  PetscTruth flg;
  PetscFunctionBegin;
  info = PetscStrcmp(str1,str2,&flg);CHKERRQ(info);
  if (flg==PETSC_TRUE) *flag=TAO_TRUE; else *flag=TAO_FALSE;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoStrCpy"
int TaoStrcpy(char* str1,const char*str2){
  int info;
  PetscFunctionBegin;
  info = PetscStrcpy(str1,str2);CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoPublish_Petsc"
static int TaoPublish_Petsc(PetscObject obj)
{
#if defined(PETSC_HAVE_AMS)
  TAO_SOLVER       v = (TAO_SOLVER) obj;
  int          info;
#endif

  TaoFunctionBegin;

#if defined(PETSC_HAVE_AMS)
  /* if it is already published then return */
  if (v->amem >=0 ) TaoFunctionReturn(0);

  info = PetscObjectPublishBaseBegin(obj);CHKERRQ(info);
  info = AMS_Memory_add_field((AMS_Memory)v->amem,"Iteration",&v->iter,1,AMS_INT,AMS_READ,AMS_COMMON,AMS_REDUCT_UNDEF);CHKERRQ(info);
  info = AMS_Memory_add_field((AMS_Memory)v->amem,"Residual",&v->fc,1,AMS_DOUBLE,AMS_READ,AMS_COMMON,AMS_REDUCT_UNDEF);CHKERRQ(info);
  info = AMS_Memory_add_field((AMS_Memory)v->amem,"Residual",&v->norm,1,AMS_DOUBLE,AMS_READ,AMS_COMMON,AMS_REDUCT_UNDEF);CHKERRQ(info);
  info = AMS_Memory_add_field((AMS_Memory)v->amem,"Iteration",(int*)&v->reason,1,AMS_INT,AMS_READ,AMS_COMMON,AMS_REDUCT_UNDEF);CHKERRQ(info);
  info = PetscObjectPublishBaseEnd(obj);CHKERRQ(info);
#endif
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoObjectCreate"
int TaoObjectCreate( TAO_SOLVER *newsolver, MPI_Comm comm){
  TAO_SOLVER solver;
  int info;

  PetscFunctionBegin;
  info = PetscHeaderCreate(solver,_p_TAO_SOLVER,int,TAO_COOKIE,-1,"TAO Solver",comm,TaoDestroy,0); CHKERRQ(info);

  ((PetscObject)solver)->bops->publish      = TaoPublish_Petsc;
  info=PetscPublishAll(solver);CHKERRQ(info);
  *newsolver = solver;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoObjectDestroy"
int TaoObjectDestroy( TAO_SOLVER solver){
  int info;
  PetscFunctionBegin;
  /* if memory was published with AMS then destroy it */

  info = PetscObjectDepublish(solver);CHKERRQ(info);
  
  PetscHeaderDestroy(solver); 

  PetscFunctionReturn(0);
}
static int one = 1;
static char *executable = (char *)"tao";
static char **executablePtr = &executable;

#undef __FUNCT__  
#define __FUNCT__ "TaoLogClassRegister"
int TaoLogClassRegister(int*cookie,const char*objectname){
  int info;
  int argc;
  char **args;
  PetscFunctionBegin;

  info = TaoGetArgs(&argc,&args); CHKERRQ(info);

#if !defined(PARCH_t3d)
  info = PetscSetHelpVersionFunctions(TaoPrintHelpIntro,TaoPrintVersion); CHKERRQ(info);
#endif

  if (!PetscInitializeCalled) {
    if (argc&&args){
      info = PetscInitialize(&argc,&args,0,0); CHKERRQ(info);
    } else {
      info = PetscInitialize(&one,&executablePtr,0,0); CHKERRQ(info);
    }
    TaoBeganPetsc = TAO_TRUE;
  }
  info=PetscCookieRegister(objectname,cookie);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#endif

/* The PetscObject is reduced and microkernal capabilities are absent
#define PetscObjectComposeFunctionDynamic(a,b,c,d) 0
#define PetscObjectQueryFunction(a,b,c) 0
  PetscLogInfo((void *,const char*,...){ return 0);}

*/
