#include "taoapp.h"  /*I  "tao.h"  I*/


#undef __FUNCT__  
#define __FUNCT__ "TaoApplicationFreeMemory"
/*@C
   TaoApplicationFreeMemory - Calls PetscFree() the argument.

   Collective on TAO_APPLICATION

   Input Parameters:
.  ctx - pointer to structure created by PetscMalloc().

   Level: intermediate

   Note:
   A pointer to this routine can be passed into TaoApp functions which will
   call this routine at a later time.

.keywords: destroy

.seealso: TaoAppSetDestroyRoutine(), PetscFree()
@*/
int TaoApplicationFreeMemory(void*ctx){
  int info; 
  PetscFunctionBegin;
  info=PetscFree(ctx); CHKERRQ(info);
  PetscFunctionReturn(0); 
}


typedef struct {
  Vec xl,xu;
} TaoVBDCtx;

#undef __FUNCT__  
#define __FUNCT__ "TaoPetscApplicationCopyBounds"
static int TaoPetscApplicationCopyBounds(TAO_APPLICATION taoapp, Vec XL, Vec XU, void*ctx){
  int info;
  TaoVBDCtx* vbctx = (TaoVBDCtx*)ctx;
  PetscFunctionBegin;
  info=VecCopy(vbctx->xl,XL);CHKERRQ(info);
  info=VecCopy(vbctx->xu,XU);CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetVariableBounds"
/*@
   TaoAppSetVariableBounds - Set bounds on the variables.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  XL - vector of lower bounds upon the solution vector
-  XU - vector of upper bounds upon the solution vector

   Note:
   This routine should be called before TaoSetApplicationSolver()

   Level: beginner

.keywords: bounds

.seealso: TaoGetVariableBoundVecs(), TaoAppSetVariableBoundsRoutine()
@*/
int TaoAppSetVariableBounds(TAO_APPLICATION taoapp, Vec XL, Vec XU){
  int info;
  TaoVBDCtx* vbctx;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(XL,VEC_COOKIE,2);
  PetscValidHeaderSpecific(XU,VEC_COOKIE,3);
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  PetscNew(TaoVBDCtx,&vbctx);
  vbctx->xl=XL;
  vbctx->xu=XU;
  info=TaoAppSetVariableBoundsRoutine(taoapp,TaoPetscApplicationCopyBounds,(void*)vbctx);
  info=TaoAppSetDestroyRoutine(taoapp,TaoApplicationFreeMemory, (void*)vbctx); CHKERRQ(info);
  PetscFunctionReturn(0);
}


