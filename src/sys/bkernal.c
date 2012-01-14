#include "tao_general.h"      /*I "tao_general.h"  I*/

#ifndef TAO_USE_PETSC

#include "src/tao_impl.h"
#include <stdlib.h>
#include <string.h>

#define MAX_SOLVER_NAME_LENGTH 20
#define MAX_NUMBER_OF_SOLVERS 20

typedef struct {
  int (*rr)(TAO_SOLVER);
  char identifier[MAX_SOLVER_NAME_LENGTH+1];
} _P_TSolver;

typedef struct {
  _P_TSolver  TSolver[MAX_NUMBER_OF_SOLVERS];
  int nsolvers;
  int maxsolvers;
  int argc;
  char **args;
} TaoSolverList;

static TaoSolverList TaoList;


#undef __FUNCT__  
#define __FUNCT__ "TaoPrintStatement"
int TaoPrintStatement(TAO_SOLVER tao, const char *statement){
  TaoFunctionBegin;
  printf(statement);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintInt"
int TaoPrintInt(TAO_SOLVER tao, const char *statement, int n){
  TaoFunctionBegin;
  printf(statement,n);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintDouble"
int TaoPrintDouble(TAO_SOLVER tao, const char *statement,double dd){
  TaoFunctionBegin;
  printf(statement,dd);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoPrintString"
int TaoPrintString(TAO_SOLVER tao, const char *statement,const char *str){
  TaoFunctionBegin;
  printf(statement,str);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionsHead"
int TaoOptionsHead(const char *heading){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionsTail"
int TaoOptionsTail(){
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionInt"
int TaoOptionInt(const char *opt,const char *text,const char *man,int defaultv,int *value,TaoTruth *set){
  TaoFunctionBegin;
  if (set) *set=TAO_FALSE;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionDouble"
int TaoOptionDouble(const char *opt,const char *text,const char *man,double defaultv,double *value,TaoTruth *set){
  TaoFunctionBegin;
  if (set) *set=TAO_FALSE;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionString"
int TaoOptionString(const char *opt,const char *text,const char *man,const char* defaultv,char *value, int len, TaoTruth *set){
  TaoFunctionBegin;
  if (set) *set=TAO_FALSE;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoOptionName"
int TaoOptionName(const char *opt,const char *text,const char *man,TaoTruth *set){
  TaoFunctionBegin;
  if (set) *set=TAO_FALSE;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoOptionList"
int TaoOptionList(const char *opt, const char *ltext, const char *man,
                  const char **list, int nlist, const char *defaultv,
                  int *value, TaoTruth *set)
{
  TaoFunctionBegin;
  if (set) *set=TAO_FALSE;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoMethodsList"
int TaoMethodsList(const char *opt,const char *ltext,const char *man,const char *defaultv,char *value,int len,TaoTruth *set){
  TaoFunctionBegin;
  if (set)  *set=TAO_FALSE;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoFindSolver"
int TaoFindSolver(TAO_SOLVER tao, TaoMethod type,  int (**r)(TAO_SOLVER) ){
  int i;
  int nsolvers=TaoList.nsolvers;
  TaoFunctionBegin;
  *r=0;
  for (i=0; i<nsolvers; i++){
    if (strncmp(type,TaoList.TSolver[i].identifier,MAX_SOLVER_NAME_LENGTH)==0){
      *r=TaoList.TSolver[i].rr;
      strncpy(tao->type_name,type,MAX_SOLVER_NAME_LENGTH);
      break;
    }
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoRegisterDestroy"
int TaoRegisterDestroy(){
  TaoFunctionBegin;
  TaoList.nsolvers=0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoRegisterDynamic"
EXTERN_C_BEGIN   // make extern c to allow extern c function pointer as argument
int TaoRegisterDynamic(const char *sname,const char *path,const char *name,int (*function)(TAO_SOLVER))
{
  int i=TaoList.nsolvers;

  TaoFunctionBegin;
  if (TaoList.nsolvers>TaoList.maxsolvers-1){
    return 1;
  } 
  TaoList.TSolver[i].rr=function;
  strncpy(TaoList.TSolver[i].identifier,sname,MAX_SOLVER_NAME_LENGTH);
  TaoList.nsolvers++;
  TaoFunctionReturn(0);
}
EXTERN_C_END

#undef __FUNCT__  
#define __FUNCT__ "TaoCompareMethod"
int TaoCompareMethod(TAO_SOLVER tao, TaoMethod type, TaoTruth *issame){
  int info;
  TaoMethod currenttype;
  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  info = TaoGetMethod(tao,&currenttype);CHKERRQ(info);
  if (currenttype && type){
    info = strncmp(type,currenttype,MAX_SOLVER_NAME_LENGTH);  
    if (info==0){
      *issame=TAO_TRUE; 
    } else {
      *issame=TAO_FALSE;
    }
  } else {
    *issame=TAO_FALSE;
  }
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoObjectCreate"
int TaoObjectCreate( TAO_SOLVER *newsolver, MPI_Comm comm){
  TAO_SOLVER solver;
  int info;
  
  TaoFunctionBegin;
  info=TaoNew(struct _p_TAO_SOLVER,(void**)&solver);CHKERRQ(info);
  //solver=(TAO_SOLVER)malloc(sizeof(struct _p_TAO_SOLVER)); 
  *newsolver = solver;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoObjectDestroy"
int TaoObjectDestroy( TAO_SOLVER solver){
  int info;
  TaoFunctionBegin;
  info=TaoFree(solver);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoStrcpy"
int TaoStrcpy(char* str1,const char*str2){
  TaoFunctionBegin;
  strcpy(str1,str2);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoStrcmp"
int TaoStrcmp(const char *str1,const char *str2,TaoTruth *flag){
  int info;
  TaoFunctionBegin;
  info = strcmp(str1,str2);
  if (info) *flag=TAO_FALSE; else *flag=TAO_TRUE;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoLogClassRegister"
int TaoLogClassRegister(int*,const char*){
  int i,info;
  TaoFunctionBegin;
  TaoList.nsolvers=0;
  TaoList.maxsolvers=MAX_NUMBER_OF_SOLVERS;
  for (i=0;i<TaoList.maxsolvers;i++){
    //    TaoList.TSolver[i].identifier[0]='';
    TaoList.TSolver[i].rr=0;
  }
  info = TaoGetArgs(&TaoList.argc,&TaoList.args);CHKERRQ(info);
  TaoFunctionReturn(0);
}

int TaoHelpPrintf(MPI_Comm,const char*,...){
  return 0;
}


#define PetscInfo  (a)   0;


#endif
