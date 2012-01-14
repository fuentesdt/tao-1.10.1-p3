/*$Id$*/

#include "private/matimpl.h"                  /*I  "mat.h"     I*/
#include "src/tao_impl.h"    /*I  "tao_solver.h"  I*/

#undef __FUNC__  
#define __FUNC__ "TaoDefaultComputeHessianColor"
/* @C
   TaoDefaultComputeHessianColor - Computes the Hessian using colored finite differences. 

   Collective on TAO_SOLVER

   Input Parameters:
+  x1 - compute Hessian at this point
-  ctx - application's gradient context, as set with TaoSetGradient()

   Output Parameters:
+  H - Hessian matrix (not altered in this routine)
.  B - newly computed Hessian matrix to use with preconditioner (generally the same as H)
-  flag - flag indicating whether the matrix sparsity structure has changed

   Options Database Keys:
.  -mat_fd_coloring_freq <freq> - Activates TaoDefaultComputeHessianColor()

   Level: advanced

 Concepts: TAO_SOLVER, finite differences, Hessian, coloring, sparse

.seealso: TaoSetHessian(), MatFDColoringCreate(), MatFDColoringSetFunction(), TaoDefaultComputeHessian()

@ */
/*
int TaoDefaultComputeHessianColor(TAO_SOLVER tao,Vec x1,Mat *H,Mat *B,MatStructure *flag,void *ctx)
{
  MatFDColoring color = (MatFDColoring) ctx;
  int           info,freq,it;

  PetscFunctionBegin;
  info = MatFDColoringGetFrequency(color,&freq);CHKERRQ(info);
  info = TaoGetIterationNumber(tao,&it);CHKERRQ(info);

  if ((freq > 1) && ((it % freq) != 1)) {
    PLogInfo(color,"TaoDefaultComputeHessianColor:Skipping Hessian, it %d, freq %d\n",it,freq);
    *flag = SAME_PRECONDITIONER;
    PetscFunctionReturn(0);
  } else {
    PLogInfo(color,"TaoDefaultComputeHessianColor:Computing Hessian, it %d, freq %d\n",it,freq);
    *flag = SAME_NONZERO_PATTERN;
  }

  info = PLogEventBegin(Tao_FunctionEval,tao,x1,0,0);
  PetscStackPush("Tao user function");
  info = MatFDColoringApply(*B,color,x1,flag,tao);CHKERRQ(info);
  PetscStackPop;
  tao->nfuncs++;
  info = PLogEventEnd(Tao_FunctionEval,tao,x1,0,0);
  PetscFunctionReturn(0);
}
*/


