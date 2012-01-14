
#include "tao_solver.h"  /*I "tao_solver.h" I*/
#include "src/tao_impl.h"

#undef __FUNCT__  
#define __FUNCT__ "TaoDestroy_UnitStep"
static int TaoDestroy_UnitStep(TAO_SOLVER tao,void *linectx)
{
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_UnitStep"
static int TaoSetOptions_UnitStep(TAO_SOLVER tao,void *linectx)
{
  int info;
  TaoFunctionBegin;
  info = TaoOptionsHead("No Unit line search options");CHKERRQ(info);
  info = TaoOptionsTail();CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoView_UnitStep"
static int TaoView_UnitStep(TAO_SOLVER tao,void *ctx)
{
  int info;
  TaoFunctionBegin;
  info=TaoPrintStatement(tao,"  Line Search: Unit Step.\n");CHKERRQ(info);
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoApply_UnitStep"
/* @ TaoApply_LineSearch - This routine takes step length of 1.0.

   Input Parameters:
+  tao - TAO_SOLVER context
.  X - current iterate (on output X contains new iterate, X + step*S)
.  S - search direction
.  f - objective function evaluated at X
.  G - gradient evaluated at X
.  W - work vector
.  gdx - inner product of gradient and the direction of the first linear manifold being searched
-  step - initial estimate of step length

   Output parameters:
+  f - objective function evaluated at new iterate, X + step*S
.  G - gradient evaluated at new iterate, X + step*S
.  X - new iterate
-  step - final step length

   Info is set to 0.

@ */
static int TaoApply_UnitStep(TAO_SOLVER tao,TaoVec* X,TaoVec* G,TaoVec* S,TaoVec* W,double *f, double *f_full,
                        double *step,TaoInt *info2,void*ctx)
{
  int       info;
  double fnew;
  TaoVec *XL,*XU;

  TaoFunctionBegin;
  tao->new_search=TAO_TRUE;
  info = W->CopyFrom(X); CHKERRQ(info);
  info = W->Axpy(*step,S);CHKERRQ(info);
  info = TaoGetVariableBounds(tao,&XL,&XU); CHKERRQ(info);
  if (XL && XU){
    info = W->Median(XL,W,XU);CHKERRQ(info);
  }
  info = TaoComputeMeritFunctionGradient(tao,W,&fnew,G); CHKERRQ(info);
  info = X->CopyFrom(W); CHKERRQ(info);
  info = PetscInfo1(tao,"Tao Apply Unit Step: %4.4e\n",*step);
         CHKERRQ(info);
  if (*f<fnew){
    info = PetscInfo2(tao,"Tao Apply Unit Step, FINCREASE: F old:= %12.10e, F new: %12.10e\n",*f,fnew); CHKERRQ(info);
  }
  *f=fnew;
  *f_full = fnew;
  *info2 = 0;
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoCreateUnitLineSearch"
/*@C
   TaoCreateUnitLineSearch - Always use step length of 1.0

   Input Parameters:
.  tao - TAO_SOLVER context

   Note:
   This routine is never used by default.

   Level: advanced

.keywords: TAO_SOLVER, linesearch
@*/
int TaoCreateUnitLineSearch(TAO_SOLVER tao)
{
  int info;

  TaoFunctionBegin;
  info = TaoSetLineSearch(tao,0,
			  TaoSetOptions_UnitStep,
			  TaoApply_UnitStep,
			  TaoView_UnitStep,
			  TaoDestroy_UnitStep,
			  0);CHKERRQ(info);

  TaoFunctionReturn(0);
}

