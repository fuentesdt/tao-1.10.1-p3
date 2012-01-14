#include "src/constrained/impls/augmented/augmented.h"

/* ---------------------------------------------------------- */
#undef __FUNC__  
#define __FUNC__ "TaoSetUp_AUGMENTED"
int TaoSetUp_AUGMENTED(TAO_SOLVER tao,void *solver)
{
  TAO_AUGMENTED *augmented = (TAO_AUGMENTED *)solver;
  int info;

  TaoFunctionBegin;
  info = TaoCheckFGH(tao);CHKERRQ(info);
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNC__  
#define __FUNC__ "TaoDestroy_AUGMENTED"
static int TaoDestroy_AUGMENTED(TAO_SOLVER tao, void *solver)
{
  TAO_AUGMENTED *augmented = (TAO_AUGMENTED *)solver;
  int info;

  TaoFunctionBegin;
  TaoFree(augmented);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNC__  
#define __FUNC__ "TaoSolve_AUGMENTED"
static int TaoSolve_AUGMENTED(TAO_SOLVER tao, void *solver)
{
  TAO_AUGMENTED *augmented = (TAO_AUGMENTED *)solver;

  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNC__  
#define __FUNC__ "TaoSetOptions_AUGMENTED"
static int TaoSetOptions_AUGMENTED(TAO_SOLVER tao, void*solver)
{
  TAO_AUGMENTED *augmented = (TAO_AUGMENTED *)solver;
  int info;

  TaoFunctionBegin;
  info = OptionsHead("Augmented Lagrangian method for constrained optimization");
         CHKERRQ(info);
  info = OptionsTail();CHKERRQ(info);
  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
#undef __FUNC__  
#define __FUNC__ "TaoView_AUGMENTED"
static int TaoView_AUGMENTED(TAO_SOLVER tao,void* solver, Viewer viewer)
{
  TAO_AUGMENTED   *augmented = (TAO_AUGMENTED *)solver;
  PetscTruth isascii;
  int info;

  TaoFunctionBegin;
  info = PetscTypeCompare((PetscObject)viewer,ASCII_VIEWER,&isascii);
         CHKERRQ(info);
  if (isascii) {
  } else {
    SETERRQ(1,1,"Viewer type not supported for this object");
  }
  TaoFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNC__  
#define __FUNC__ "TaoCreate_AUGMENTED"
int TaoCreate_AUGMENTED(TAO_SOLVER tao)
{
  TAO_AUGMENTED *augmented;
  int info;

  TaoFunctionBegin;

  augmented = TaoNew(TAO_AUGMENTED); CHKPTRQ(augmented);
  PLogObjectMemory(tao, sizeof(TAO_AUGMENTED));

  info = TaoSetSolver(tao,TaoSetUp_AUGMENTED,TaoSetOptions_AUGMENTED,
                      TaoSolve_AUGMENTED,TaoView_AUGMENTED,
		      TaoDestroy_AUGMENTED,(void*)augmented); CHKERRQ(info);

  TaoFunctionReturn(0);
}
EXTERN_C_END
