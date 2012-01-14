
/*
  Implements the default TAO approach for computing the differencing
  parameter, h, used with the finite difference based matrix-free 
  Hessian-vector products.

  To make your own: clone this file and modify for your needs.

  Mandatory functions:
  -------------------
      MatTaoMFCompute_  - for a given point and direction computes h

      MatTaoMFCreate_ - fills in the MatTaoMF data structure
                           for this particular implementation

      
   Optional functions:
   -------------------
      MatTaoMFView_ - prints information about the parameters being used.
                       This is called when TaoView() or -tao_view is used.

      MatTaoMFPrintHelp_ - prints a help message on what options are
                          available for this implementation

      MatTaoMFSetFromOptions_ - checks the options database for options that 
                               apply to this method.

      MatTaoMFDestroy_ - frees any space allocated by the routines above

*/

/*
    This include file defines the data structure  MatTaoMF that 
   includes information about the computation of h. It is shared by 
   all implementations that people provide
*/
#include "src/supplementary/matrix/matrixfree/tao_mfj.h"   /*I  "tao_solver.h"   I*/

/*
   The default method has one parameter that is used to 
   "cutoff" very small values. This is stored in a data structure
   that is only visible to this file. If your method has no parameters
   it can omit this, if it has several simply reorganize the data structure.
   The data structure is "hung-off" the MatTaoMF data structure in
   the void *hctx; field.
*/
typedef struct {
  double umin;          /* minimum allowable u'a value relative to |u|_1 */
} MatTaoMFDefault;

#undef __FUNC__  
#define __FUNC__ "MatTaoMFCompute_Default"
/*
   MatTaoMFCompute_Default - Code for computing the differencing 
   paramter (h) for use with matrix-free finite differences.

   Input Parameters:
+  ctx - the matrix free context
.  U - the location at which you want the Hessian
-  a - the direction you want the derivative

   Output Parameter:
.  h - the scale computed

*/
static int MatTaoMFCompute_Default(MatTaoMFCtx ctx,Vec U,Vec a,Scalar *h)
{
  MatTaoMFDefault *hctx = (MatTaoMFDefault *) ctx->hctx;
  double           norm, sum, umin = hctx->umin;
  Scalar           dot;
  int              info;

  PetscFunctionBegin;
  /*
     This algorithm requires 2 norms and 1 inner product. Rather than
     use directly the VecNorm() and VecDot() routines (and thus have 
     three separate collective operations, we use the VecxxxBegin/End() routines
     and manually call MPI for the collective phase.

  */
  info = VecDotBegin(U,a,&dot);CHKERRQ(info);
  info = VecNormBegin(a,NORM_1,&sum);CHKERRQ(info);
  info = VecNormBegin(a,NORM_2,&norm);CHKERRQ(info);
  info = VecDotEnd(U,a,&dot);CHKERRQ(info);
  info = VecNormEnd(a,NORM_1,&sum);CHKERRQ(info);
  info = VecNormEnd(a,NORM_2,&norm);CHKERRQ(info);

  /* 
     Safeguard for step sizes that are "too small"
  */
  if (sum == 0.0) {dot = 1.0; norm = 1.0;}
#if defined(PETSC_USE_COMPLEX)
  else if (PetscAbsScalar(dot) < umin*sum && PetscReal(dot) >= 0.0) dot = umin*sum;
  else if (PetscAbsScalar(dot) < 0.0 && PetscReal(dot) > -umin*sum) dot = -umin*sum;
#else
  else if (dot < umin*sum && dot >= 0.0) dot = umin*sum;
  else if (dot < 0.0 && dot > -umin*sum) dot = -umin*sum;
#endif
  *h = ctx->error_rel*dot/(norm*norm);
  PetscFunctionReturn(0);
} 

#undef __FUNC__  
#define __FUNC__ "MatTaoMFView_Default"
/*
   MatTaoMFView_Default - Prints information about this particular 
   method for computing h. Note that this does not print the general
   information about the matrix-free method, as such info is printed
   by the calling routine.

   Input Parameters:
+  ctx - the matrix free context
-  viewer - the PETSc viewer
*/   
static int MatTaoMFView_Default(MatTaoMFCtx ctx,Viewer viewer)
{
  FILE            *fd;
  MatTaoMFDefault *hctx = (MatTaoMFDefault *)ctx->hctx;
  int             info;
  TaoTruth      isascii;

  PetscFunctionBegin;
  info = ViewerASCIIGetPointer(viewer,&fd);CHKERRQ(info);
  
  /*
     Currently this only handles the ascii file viewers, others
     could be added, but for this type of object other viewers
     make less sense
  */
  info = PetscTypeCompare((PetscObject)viewer,ASCII_VIEWER,&isascii);CHKERRQ(info);
  if (isascii) {
    info = PetscFPrintf(ctx->comm,fd,"    umin=%g (minimum iterate parameter)\n",hctx->umin);CHKERRQ(info); 
  } else {
    SETERRQ(1,1,"Viewer type not supported for this object");
  }    
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFPrintHelp_Default"
/*
   MatTaoMFPrintHelp_Default - Prints a list of all the options 
   this particular method supports.

   Input Parameter:
.  ctx - the matrix free context

*/
static int MatTaoMFPrintHelp_Default(MatTaoMFCtx ctx)
{
  char*            p;
  MatTaoMFDefault *hctx = (MatTaoMFDefault *)ctx->hctx;
  int              info;

  PetscFunctionBegin;
  info = PetscObjectGetOptionsPrefix((PetscObject)ctx->tao,&p);CHKERRQ(info);
  if (!p) p = "";
  info = (*PetscHelpPrintf)(ctx->comm,"   -%stao_mf_umin <umin> see users manual (default %g)\n",p,hctx->umin);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFSetFromOptions_Default"
/*
   MatTaoMFSetFromOptions_Default - Looks in the options database for 
   any options appropriate for this method.

   Input Parameter:
.  ctx - the matrix free context

*/
static int MatTaoMFSetFromOptions_Default(MatTaoMFCtx ctx)
{
  char*      p;
  int        info;
  double     umin;
  TaoTruth flg;

  PetscFunctionBegin;
  info = PetscObjectGetOptionsPrefix((PetscObject)ctx->tao,&p);CHKERRQ(info);
  info = OptionsGetDouble(p,"-tao_mf_umin",&umin,&flg);CHKERRQ(info);
  if (flg) {
    info = MatTaoMFDefaultSetUmin(ctx->mat,umin);CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFDestroy_Default"
/*
   MatTaoMFDestroy_Default - Frees the space allocated by 
   MatTaoMFCreate_Default(). 

   Input Parameter:
.  ctx - the matrix free context

   Notes: 
   Does not free the ctx, that is handled by the calling routine
*/
static int MatTaoMFDestroy_Default(MatTaoMFCtx ctx)
{
  int info;

  PetscFunctionBegin;
  info = PetscFree(ctx->hctx);CHKERRQ(info);
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNC__  
#define __FUNC__ "MatTaoMFDefaultSetUmin_Private"
/*
   The following two routines use the PetscObjectCompose() and PetscObjectQuery()
   mechanism to allow the user to change the Umin parameter used in this method.
*/
int MatTaoMFDefaultSetUmin_Private(Mat mat,double umin)
{
  MatTaoMFCtx     ctx;
  MatTaoMFDefault *hctx;
  int              info;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  if (!ctx) {
    SETERRQ(1,1,"MatTaoMFDefaultSetUmin() attached to non-shell matrix");
  }
  hctx = (MatTaoMFDefault *) ctx->hctx;
  hctx->umin = umin;

 PetscFunctionReturn(0);
} 
EXTERN_C_END

#undef __FUNC__  
#define __FUNC__ "MatTaoMFDefaultSetUmin"
/*@
    MatTaoMFDefaultSetUmin - Sets the "umin" parameter used by the default
    TAO routine for computing the differencing parameter, h, which is used
    for matrix-free Hessian-vector products that are based on the gradient.

   Input Parameters:
+  A - the matrix created with MatCreateTaoMF()
-  umin - the parameter

   Level: advanced

   Notes:
   See the manual page for MatCreateTaoMF() for a complete description of the
   algorithm used to compute h.

.seealso: MatTaoMFSetFunctionError(), MatCreateTaoMF()
@*/
int MatTaoMFDefaultSetUmin(Mat A,double umin)
{
  int info, (*f)(Mat,double);

  PetscFunctionBegin;
  PetscValidHeaderSpecific(A,MAT_COOKIE);
  info = PetscObjectQueryFunction((PetscObject)A,"MatTaoMFDefaultSetUmin_C",(void **)&f);CHKERRQ(info);
  if (f) {
    info = (*f)(A,umin);CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNC__  
#define __FUNC__ "MatTaoMFCreate_Default"
/*
   MatTaoMFCreate_Default - Standard PETSc code for 
   computing h with matrix-free finite differences.

   Input Parameter:
.  ctx - the matrix free context created by MatTaoMFCreate()

*/
int MatTaoMFCreate_Default(MatTaoMFCtx ctx)
{
  MatTaoMFDefault *hctx;
  int              info;

  PetscFunctionBegin;

  /* allocate my own private data structure */
  hctx                     = (MatTaoMFDefault *)PetscMalloc(sizeof(MatTaoMFDefault));CHKPTRQ(hctx);
  ctx->hctx                = (void *) hctx;
  /* set a default for my parameter */
  hctx->umin               = 1.e-6;

  /* set the functions I am providing */
  ctx->ops->compute        = MatTaoMFCompute_Default;
  ctx->ops->destroy        = MatTaoMFDestroy_Default;
  ctx->ops->view           = MatTaoMFView_Default;  
  ctx->ops->printhelp      = MatTaoMFPrintHelp_Default;  
  ctx->ops->setfromoptions = MatTaoMFSetFromOptions_Default;  

  info = PetscObjectComposeFunctionDynamic((PetscObject)ctx->mat,"MatTaoMFDefaultSetUmin_C",
                            "MatTaoMFDefaultSetUmin_Private",
                            (void *) MatTaoMFDefaultSetUmin_Private);CHKERRQ(info);
  PetscFunctionReturn(0);
}
EXTERN_C_END







