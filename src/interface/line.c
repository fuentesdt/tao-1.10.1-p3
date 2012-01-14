/*$Id$*/

#include "src/tao_impl.h"      /*I "tao_solver.h"  I*/

#undef __FUNCT__  
#define __FUNCT__ "TaoSetLineSearch"
/*@C
   TaoSetLineSearch - Set the line search routine for algorithms that 
   require one.

   Collective on TAO_SOLVER

   Input Parameter:
+  tao - the TAO_SOLVER solver context
.  setup -  setup routine (or TAO_NULL)
.  options - set line search options (or TAO_NULL)
.  line - the line search routine
.  viewit - routine that views the linesearch (or TAO_NULL)
.  destroy - destroys the user defined routine when the solver is destroyed (or TAO_NULL)
-  ctx - linesearch structure (or TAO_NULL)

   Calling sequence of line:
$     line(TAO_SOLVER tao, TaoVec* X, TaoVec* G, TaoVec* DX, TaoVec* Work, double *f,double *step, int *flag, void *ctx)

   Input Parameter for line search:
+  tao - the TAO_SOLVER solver context
.  xx - current solution
.  gg - current gradient
.  dxdx - step direction
.  Work - work vector
.  f - function value
.  step - initial stepsize
-  ctx - user-defined line search context 

   Output Parameter for line search:
+  X - new solution
.  G - new gradient
.  f - new function value
.  step - multiple of DX added to the previous solution
-  flag - indicator of success or failure (flag=0 is a success, flag=7 means DX is not a descent direction)

   Notes:
   The input parameter gdx should be negative and is used to test the Armijo 
   condition.
   
   To ensure descent in a projected line search, the input parameter gdx 
   should be the inner product of the gradient and the first linear manifold 
   being searched.

   Level: advanced

.keywords: TAO_SOLVER, destroy

.seealso: TaoLineSearchApply()
@*/
int TaoSetLineSearch(TAO_SOLVER tao, 
		     int (*setup)(TAO_SOLVER,void*),
		     int (*options)(TAO_SOLVER,void*),
		     int (*line)(TAO_SOLVER,TaoVec*,TaoVec*,TaoVec*,TaoVec*,
				 double*,double*,double*,TaoInt*,void*),
		     int (*viewit)(TAO_SOLVER,void*),
		     int (*destroy)(TAO_SOLVER,void*),
		     void *ctx)
{
  int info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);

  info=TaoLineSearchDestroy(tao);CHKERRQ(info);

  tao->LineSearchSetUp	        = setup;
  tao->LineSearchSetFromOptions = options; 
  tao->LineSearchApply          = line; 
  tao->LineSearchView           = viewit;
  tao->LineSearchDestroy        = destroy;
  tao->linectx                  = ctx;
  tao->lsflag                  = 0;

  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoLineSearchView"
/*@C
   TaoLineSearchView - View information about the line search

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Level: developer

.keywords: line search, view

.seealso: TaoView()
@*/
int TaoLineSearchView(TAO_SOLVER tao)
{
  int info;
  
  if (tao->LineSearchView){
    info = (*tao->LineSearchView)(tao,tao->linectx);CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoLineSearchSetUp"
/*@C
   TaoLineSearchSetUp - Setup the line search for an optimization solver

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Level: developer

.keywords: line search, setup

.seealso: TaoSetUp()
@*/
int TaoLineSearchSetUp(TAO_SOLVER tao)
{
  int info;

  if (tao->LineSearchSetUp) {
    info = (*tao->LineSearchSetUp)(tao,tao->linectx);CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoLineSearchDestroy"
/*@C
   TaoLineSearchDestroy - Destroy the line search in an TAO solver

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Level: developer

.keywords: line search, destroy

.seealso: TaoDestroy()
@*/
int TaoLineSearchDestroy(TAO_SOLVER tao)
{
  int info;

  if (tao->LineSearchDestroy && tao->linectx){
    info = (*tao->LineSearchDestroy)(tao,tao->linectx);CHKERRQ(info);
  }
  tao->linectx=0;
  TaoFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoLineSearchApply"
/*@C
   TaoLineSearchApply - Applies a line search.

   Collective on TAO_SOLVER

   Input Parameter:
+  tao - the TAO_SOLVER solver context
.  xx - current solution
.  gg - current gradient
.  dxdx - step direction
.  ww - work vector
.  f - function value
-  step - initial stepsize

   Output Parameter:
+  xx - new solution
.  gg - new gradient
.  f - new function value
.  step - multiple of DX added to the previous solution
-  flag - indicator of success or failure (flag=0 is a success, flag=7 means DX is not a descent direction)

   Notes:
   The input parameter gdx should be negative and is used to test the Armijo 
   condition.
   
   To ensure descent in a projected line search, the input parameter gdx 
   should be the inner product of the gradient and the first linear manifold 
   being searched.

   Level: developer

.keywords: line search, help

.seealso: TaoSetLineSearch()
@*/
int TaoLineSearchApply(TAO_SOLVER tao, TaoVec *xx, TaoVec *gg, TaoVec *dxdx, TaoVec *ww, 
		       double *f, double *f_full, double *step, TaoInt*flag)
{
  int info;

  TaoFunctionBegin;
  TaoValidHeaderSpecific(tao,TAO_COOKIE,1);
  if (tao->LineSearchApply){
    info = (*tao->LineSearchApply)(tao,xx,gg,dxdx,ww,f,f_full,step,flag,tao->linectx);
    CHKERRQ(info);
    tao->lsflag=*flag;
  } else {
    SETERRQ(1,"No line search has been set.");
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoLineSearchSetFromOptions"
/*@C
   TaoLineSearchSetFromOptions - Set options for the line search in an optimization solver

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Level: developer

.keywords: line search, options

.seealso: TaoSetFromOptions()
@*/
int TaoLineSearchSetFromOptions(TAO_SOLVER tao)
{
  int info;

  if (tao->LineSearchSetFromOptions){
    info = (*tao->LineSearchSetFromOptions)(tao,tao->linectx); CHKERRQ(info);
  }
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoGetCurrentStepLength"
/*@
   TaoGetCurrentStepLength -- Get the current step length.  Only valid during line searches.

   Collective on TAO_SOLVER

   Input Parameter:
.  tao - the TAO_SOLVER solver context

   Output Parameter:
+  steplength - Current steplength candidate in line search
-  newsearch - TAO_TRUE if this is the start of a new search, TAO_FALSE otherwise.  Use TAO_NULL to ignore.

   Level: advanced

.keywords: line search,
@*/
int TaoGetCurrentStepLength(TAO_SOLVER tao, double *steplength, TaoTruth *newsearch)
{
  *steplength = tao->current_step;
  if (newsearch) *newsearch = tao->new_search;
  return 0;
}

