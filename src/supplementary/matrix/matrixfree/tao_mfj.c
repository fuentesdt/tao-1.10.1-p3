/*$Id$*/

#include "src/tao_impl.h"
#include "src/supplementary/matrix/matrixfree/tao_mfj.h" /*I  "tao_solver.h"   I*/

FList MatTaoMFList              = 0;
int   MatTaoMFRegisterAllCalled = 0;

#undef __FUNC__  
#define __FUNC__ "MatTaoMFSetType"
/*@C
    MatTaoMFSetType - Sets the method that is used to compute the 
    differencing parameter for finite difference matrix-free formulations. 

    Input Parameters:
+   mat - the "matrix-free" matrix created via MatCreateTaoMF()
-   ftype - the type requested

    Level: advanced

    Notes:
    For example, such routines can compute h for use in
    Hessian-vector products of the form

                      g(u+ha) - g(u)
          H(u)a  ~=  ---------------- ,
                              h

    where g(u) denotes the gradient of f(u), and H(u) denotes
    the Hessian of f(u).

.seealso: MatCreateTaoMF(), MatTaoMFRegister()
@*/
int MatTaoMFSetType(Mat mat,char *ftype)
{
  int         info, (*r)(MatTaoMFCtx);
  MatTaoMFCtx ctx;
  TaoTruth  match;
  
  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);

  /* already set, so just return */
  info = PetscTypeCompare((PetscObject)ctx,ftype,&match);CHKERRQ(info);
  if (match) PetscFunctionReturn(0);


  /* destroy the old one if it exists */
  if (ctx->ops->destroy) {
    info = (*ctx->ops->destroy)(ctx);CHKERRQ(info);
  }

  /* Get the function pointers for the requrested method */
  if (!MatTaoMFRegisterAllCalled) {
    info = MatTaoMFRegisterAll("${TAO_DIR}/lib/${PETSC_ARCH}/libtao.so");CHKERRQ(info);
  }

  info =  FListFind(ctx->comm, MatTaoMFList, ftype,(int (**)(void *)) &r );CHKERRQ(info);

  if (!r) SETERRQ(1,1,"Unknown MatTaoMF type given");

  info = (*r)(ctx);CHKERRQ(info);
  info = PetscStrncpy(ctx->type_name,ftype,256);CHKERRQ(info);

  PetscFunctionReturn(0);
}

/*MC
   MatTaoMFRegisterDynamic - Adds a method to the MatTaoMF registry.

   Synopsis:
   MatTaoMFRegisterDynamic(char *name_solver,char *path,char *name_create,int (*routine_create)(MatTaoMF))

   Not Collective

   Input Parameters:
+  name_solver - name of a new user-defined compute-h module
.  path - path (either absolute or relative) the library containing this solver
.  name_create - name of routine to create method context
-  routine_create - routine to create method context

   Level: developer

   Notes:
   MatTaoMFRegisterDynamic() may be called multiple times to add several user-defined solvers.

   If dynamic libraries are used, then the fourth input argument (routine_create)
   is ignored.

   Environmental variables such as ${TAO_DIR}, ${PETSC_ARCH}, ${PETSC_DIR}
   and others of the form ${any_environmental_variable} occuring in pathname will be 
   replaced with appropriate values.

   Sample usage:
.vb
   MatTaoMFRegisterDynamic("my_h",/home/username/my_lib/lib/solaris/mylib.a,
               "MyHCreate",MyHCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     MatTaoMFSetType(mfctx,"my_h")
   or at runtime via the option
$     -tao_mf_type my_h

.keywords: MatTaoMF, register

.seealso: MatTaoMFRegisterAll(), MatTaoMFRegisterDestroy()
M*/

#undef __FUNC__  
#define __FUNC__ "MatTaoMFRegister"
int MatTaoMFRegister(char *sname,char *path,char *name,int (*function)(MatTaoMFCtx))
{
  int info;
  char fullname[256];

  PetscFunctionBegin;
  info = FListConcat(path,name,fullname); CHKERRQ(info);
  info = FListAdd(&MatTaoMFList,sname,fullname,(int (*)(void*))function);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFRegisterDestroy"
/*@C
   MatTaoMFRegisterDestroy - Frees the list of MatTaoMF methods that were
   registered by MatTaoMFRegister().

   Not Collective

   Level: developer

.keywords: MatTaoMF, register, destroy

.seealso: MatTaoMFRegister(), MatTaoMFRegisterAll()
@*/
int MatTaoMFRegisterDestroy(void)
{
  int info;

  PetscFunctionBegin;
  if (MatTaoMFList) {
    info = FListDestroy( &MatTaoMFList );CHKERRQ(info);
    /*    info = FListDestroy( MatTaoMFList );CHKERRQ(info); */
    MatTaoMFList = 0;
  }
  MatTaoMFRegisterAllCalled = 0;
  PetscFunctionReturn(0);
}

/* ----------------------------------------------------------------------------------------*/
#undef __FUNC__  
#define __FUNC__ "MatTaoMFDestroy_Private"
int MatTaoMFDestroy_Private(Mat mat)
{
  int          info;
  MatTaoMFCtx ctx;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  info = VecDestroy(ctx->w);CHKERRQ(info);
  if (ctx->ops->destroy) {info = (*ctx->ops->destroy)(ctx);CHKERRQ(info);}
  if (ctx->sp) {info = MatNullSpaceDestroy(ctx->sp);CHKERRQ(info);}
  info = PetscFree(ctx->ops);CHKERRQ(info);
  info = PetscFree(ctx);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFView_Private"
/*
   MatTaoMFView_Private - Views matrix-free parameters.

*/
int MatTaoMFView_Private(Mat J,Viewer viewer)
{
  int         info;
  MatTaoMFCtx ctx;
  MPI_Comm    comm;
  FILE        *fd=0;
  TaoTruth  isascii;

  PetscFunctionBegin;
  info = PetscObjectGetComm((PetscObject)J,&comm);CHKERRQ(info);
  info = MatShellGetContext(J,(void **)&ctx);CHKERRQ(info);
  info = PetscTypeCompare((PetscObject)viewer,ASCII_VIEWER,&isascii);CHKERRQ(info);
  if (isascii) {
     info = PetscFPrintf(comm,fd,"  TAO_SOLVER matrix-free approximation:\n");CHKERRQ(info);
     info = PetscFPrintf(comm,fd,"    err=%g (relative error in function evaluation)\n",ctx->error_rel);CHKERRQ(info);
     info = PetscFPrintf(ctx->comm,fd,"    Using %s compute h routine\n",ctx->type_name);CHKERRQ(info);
     if (ctx->ops->view) {
       info = (*ctx->ops->view)(ctx,viewer);CHKERRQ(info);
     }
  } else {
    SETERRQ(1,1,"Viewer type not supported for this object");
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFAssemblyEnd_Private"
/*
   MatTaoMFAssemblyEnd_Private - Resets the ctx->ncurrenth to zero. This 
   allows the user to indicate the beginning of a new linear solve by calling
   MatAssemblyXXX() on the matrix free matrix. This then allows the 
   MatTaoMFCreate_WP() to properly compute ||U|| only the first time
   in the linear solver rather than every time.
*/
int MatTaoMFAssemblyEnd_Private(Mat J)
{
  int            info;

  PetscFunctionBegin;
  info = MatTaoMFResetHHistory(J);CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNC__  
#define __FUNC__ "MatTaoMFMult_Private"
/*
  MatTaoMFMult_Private - Default matrix-free form for Hessian-vector
  product, y = H(u)*a, where H(u) denotes the Hessian of f(u).

        y ~= ( g(u + ha) - g(u) )/h, 
  where g(u) = gradient of f(u), as set by TaoSetGradient()
        u    = current iterate
        h    = difference interval
*/
int MatTaoMFMult_Private(Mat mat,Vec a,Vec y)
{
  MatTaoMFCtx ctx;
  TAO_SOLVER  tao;
  Scalar      h, mone = -1.0;
  Vec         w, U, F;
  int         info, (*eval_fct)(TAO_SOLVER,Vec,Vec)=0;

  PetscFunctionBegin;
  /* We log matrix-free matrix-vector products separately, so that we can
     separate the performance monitoring from the cases that use conventional
     storage.  We may eventually modify event logging to associate events
     with particular objects, hence alleviating the more general problem. */
  info = PLogEventBegin(MAT_MatrixFreeMult,a,y,0,0);

  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  tao  = ctx->tao;
  w    = ctx->w;
  info = TaoGetSolution(tao,&U);CHKERRQ(info);

  /* 
      Compute differencing parameter 
  */
  if (!ctx->ops->compute) {
    info = MatTaoMFSetType(mat,"default");CHKERRQ(info);
    info = MatTaoMFSetFromOptions(mat);CHKERRQ(info);
  }
  info = (*ctx->ops->compute)(ctx,U,a,&h);CHKERRQ(info);

  /* keep a record of the current differencing parameter h */  
  ctx->currenth = h;
#if defined(PETSC_USE_COMPLEX)
  PLogInfo(mat,"Current differencing parameter: %g + %g i\n",PetscReal(h),PetscImaginary(h));
#else
  PLogInfo(mat,"Current differencing parameter: %g\n",h);
#endif
  if (ctx->historyh && ctx->ncurrenth < ctx->maxcurrenth) {
    ctx->historyh[ctx->ncurrenth] = h;
  }
  ctx->ncurrenth++;

  /* w = u + ha */
  info = VecWAXPY(w,h,a,U);CHKERRQ(info);


  if (!ctx->func) {
    eval_fct = TaoComputeGradient;
    info     = TaoGetGradient(tao,&F);CHKERRQ(info);
    info = eval_fct(tao,w,y);CHKERRQ(info);
  } else {
    F = ctx->funcvec;
    /* compute func(U) as base for differencing */
    if (ctx->ncurrenth == 1) {
      info = (*ctx->func)(tao,U,F,ctx->funcctx);CHKERRQ(info);
    }
    info = (*ctx->func)(tao,w,y,ctx->funcctx);CHKERRQ(info);
  }

  info = VecAXPY(y, mone, F);CHKERRQ(info);
  h    = 1.0/h;
  info = VecScale(y, h);CHKERRQ(info);
  if (ctx->sp) {info = MatNullSpaceRemove(ctx->sp,y,PETSC_NULL);CHKERRQ(info);}

  info = PLogEventEnd(MAT_MatrixFreeMult,a,y,0,0);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatCreateTaoMF"
/*@C
   MatCreateTaoMF - Creates a matrix-free matrix context for use with
   a TAO_SOLVER solver.  This matrix can be used as the Hessian argument for
   the routine TaoSetHessian().

   Collective on TAO_SOLVER and Vec

   Input Parameters:
+  tao - the TAO_SOLVER context
-  x - vector where solution is to be stored.

   Output Parameter:
.  J - the matrix-free matrix

   Level: advanced

   Notes:
   The matrix-free matrix context merely contains the function pointers
   and work space for performing finite difference approximations of
   Hessian-vector products, H(u)*a, where H(u) denotes the Hessian
   of f(u).

   The default code uses the following approach to compute h

.vb
     H(u)*a = [g(u+h*a) - g(u)] / h,
     where
     g(u) denotes the gradient of f(u),
     H(u) denotes the Hessian of f(u), and

     h = error_rel*u'a/||a||^2                        if  |u'a| > umin*||a||_{1}
       = error_rel*umin*sign(u'a)*||a||_{1}/||a||^2   otherwise
 where
     error_rel = square root of relative error in gradient evaluation
     umin = minimum iterate parameter
.ve

   The user can set the error_rel via MatTaoMFSetGradientError() and 
   umin via MatTaoMFDefaultSetUmin().

   The user should call MatDestroy() when finished with the matrix-free
   matrix context.

   Options Database Keys:
+  -tao_mf_err <error_rel> - Sets error_rel
.  -tao_mf_unim <umin> - Sets umin (for default PETSc routine that computes h only)
-  -tao_mf_ksp_monitor - KSP monitor routine that prints differencing h

.keywords: TAO_SOLVER, default, matrix-free, create, matrix

.seealso: MatDestroy(), MatTaoMFSetGradientError(), MatTaoMFDefaultSetUmin()
          MatTaoMFSetHHistory(), MatTaoMFResetHHistory(),
          MatTaoMFGetH(),MatTaoMFKSPMonitor(), MatTaoMFRegister()
 
@*/
int MatCreateTaoMF(TAO_SOLVER tao,Vec x, Mat *J)
{
  MPI_Comm     comm;
  MatTaoMFCtx mfctx;
  int          n, nloc, info;

  PetscFunctionBegin;
  mfctx = (MatTaoMFCtx) PetscMalloc(sizeof(struct _p_MatTaoMFCtx));CHKPTRQ(mfctx);
  PLogObjectMemory(tao,sizeof(MatTaoMFCtx));
  mfctx->comm         = tao->comm;
  mfctx->sp           = 0;
  mfctx->tao          = tao;
  mfctx->error_rel    = 1.e-8; /* assumes double precision */
  mfctx->currenth     = 0.0;
  mfctx->historyh     = TAO_NULL;
  mfctx->ncurrenth    = 0;
  mfctx->maxcurrenth  = 0;
  info = PetscMemzero(mfctx->type_name,256*sizeof(char));CHKERRQ(info);

  /* 
     Create the empty data structure to contain compute-h routines.
     These will be filled in below from the command line options or 
     a later call with MatTaoMFSetType() or if that is not called 
     then it will default in the first use of MatTaoMFMult_private()
  */
  mfctx->ops                 = (TaoMFOps *)PetscMalloc(sizeof(TaoMFOps));CHKPTRQ(mfctx->ops); 
  mfctx->ops->compute        = 0;
  mfctx->ops->destroy        = 0;
  mfctx->ops->view           = 0;
  mfctx->ops->printhelp      = 0;
  mfctx->ops->setfromoptions = 0;
  mfctx->hctx                = 0;

  mfctx->func                = 0;
  mfctx->funcctx             = 0;
  mfctx->funcvec             = 0;

  info = VecDuplicate(x,&mfctx->w);CHKERRQ(info);
  info = PetscObjectGetComm((PetscObject)x,&comm);CHKERRQ(info);
  info = VecGetSize(x,&n);CHKERRQ(info);
  info = VecGetLocalSize(x,&nloc);CHKERRQ(info);
  info = MatCreateShell(comm,nloc,nloc,n,n,mfctx,J);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_MULT,(void*)MatTaoMFMult_Private);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_DESTROY,(void *)MatTaoMFDestroy_Private);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_VIEW,(void *)MatTaoMFView_Private);CHKERRQ(info);
  info = MatShellSetOperation(*J,MATOP_ASSEMBLY_END,(void *)MatTaoMFAssemblyEnd_Private);CHKERRQ(info);
  PLogObjectParent(*J,mfctx->w);
  PLogObjectParent(tao,*J);

  mfctx->mat = *J;

  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFSetFromOptions"
/*@
   MatTaoMFSetFromOptions - Sets the MatTaoMF options from the command line
   parameter.

   Collective on Mat

   Input Parameters:
.  mat - the matrix obtained with MatCreateTaoMF()

   Options Database Keys:
+  -tao_mf_type - <default,wp>
-  -tao_mf_err - square root of estimated relative error in function evaluation

   Level: advanced

.keywords: TAO_SOLVER, matrix-free, parameters

.seealso: MatCreateTaoMF(),MatTaoMFSetHHistory(), 
          MatTaoMFResetHHistory(), MatTaoMFKSPMonitor()
@*/
int MatTaoMFSetFromOptions(Mat mat)
{
  MatTaoMFCtx mfctx;
  int         info;
  char        ftype[256],p[64];
  TaoTruth  flg;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&mfctx);CHKERRQ(info);
  if (mfctx) {
    /* allow user to set the type */
    info = OptionsGetString(mfctx->tao->prefix,"-tao_mf_type",ftype,256,&flg);CHKERRQ(info);
    if (flg) {
      info = MatTaoMFSetType(mat,ftype);CHKERRQ(info);
    }

    info = OptionsGetDouble(mfctx->tao->prefix,"-tao_mf_err",&mfctx->error_rel,&flg);CHKERRQ(info);
    if (mfctx->ops->setfromoptions) {
      info = (*mfctx->ops->setfromoptions)(mfctx);CHKERRQ(info);
    }

    info = OptionsHasName(TAO_NULL,"-help",&flg);CHKERRQ(info);
    info = PetscStrcpy(p,"-");CHKERRQ(info);
    if (mfctx->tao->prefix) {info = PetscStrcat(p,mfctx->tao->prefix);CHKERRQ(info);}
    if (flg) {
      info = (*PetscHelpPrintf)(mfctx->tao->comm,"   %stao_mf_err <err>: set sqrt rel error in function (default %g)\n",p,mfctx->error_rel);CHKERRQ(info);
      if (mfctx->ops->printhelp) {
        info = (*mfctx->ops->printhelp)(mfctx);CHKERRQ(info);
      }
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFGetH"
/*@
   MatTaoMFGetH - Gets the last value that was used as the differencing 
   parameter.

   Not Collective

   Input Parameters:
.  mat - the matrix obtained with MatCreateTaoMF()

   Output Paramter:
.  h - the differencing step size

   Level: advanced

.keywords: TAO_SOLVER, matrix-free, parameters

.seealso: MatCreateTaoMF(),MatTaoMFSetHHistory(), 
          MatTaoMFResetHHistory(),MatTaoMFKSPMonitor()
@*/
int MatTaoMFGetH(Mat mat,Scalar *h)
{
  MatTaoMFCtx ctx;
  int          info;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  if (ctx) {
    *h = ctx->currenth;
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFKSPMonitor"
/*
   MatTaoMFKSPMonitor - A KSP monitor for use with the default PETSc
   TAO_SOLVER matrix free routines. Prints the differencing parameter used at 
   each step.
*/
int MatTaoMFKSPMonitor(KSP ksp,int n,double rnorm,void *dummy)
{
  PC             pc;
  MatTaoMFCtx   ctx;
  int            info;
  Mat            mat;
  MPI_Comm       comm;
  TaoTruth     nonzeroinitialguess;

  PetscFunctionBegin;
  info = PetscObjectGetComm((PetscObject)ksp,&comm);CHKERRQ(info);
  info = KSPGetPC(ksp,&pc);CHKERRQ(info);
  info = KSPGetInitialGuessNonzero(ksp,&nonzeroinitialguess);CHKERRQ(info);
  info = PCGetOperators(pc,&mat,TAO_NULL,TAO_NULL);CHKERRQ(info);
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  if (!ctx) {
    SETERRQ(1,1,"Matrix is not a matrix free shell matrix");
  }
  if (n > 0 || nonzeroinitialguess) {
#if defined(PETSC_USE_COMPLEX)
    info = PetscPrintf(comm,"%d KSP Residual norm %14.12e h %g + %g i\n",n,rnorm,
                PetscReal(ctx->currenth),PetscImaginary(ctx->currenth));CHKERRQ(info);
#else
    info = PetscPrintf(comm,"%d KSP Residual norm %14.12e h %g \n",n,rnorm,ctx->currenth);CHKERRQ(info); 
#endif
  } else {
    info = PetscPrintf(comm,"%d KSP Residual norm %14.12e\n",n,rnorm);CHKERRQ(info); 
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFSetGradient"
/*@C
   MatTaoMFSetGradient - Sets the gradient evaluation routine for use 
   in applying matrix-free Hessian-vector products.

   Collective on Mat

   Input Parameters:
+  mat     - the matrix-free matrix created via MatCreateTaoMF()
.  v       - workspace vector
.  func    - the function to use for gradient evaluation
-  funcctx - optional gradient-evaluation context passed to function

   Level: advanced

   Notes:
   If you call MatTaoMFSetGradient(), then you MUST call MatAssemblyBegin()/MatAssemblyEnd() 
   on the matrix-free matrix inside your Hessian computation routine.

   If this routine is not called, then the code will use the function that was
   set with TaoSetGradient().

.keywords: TAO_SOLVER, matrix-free, gradient

.seealso: MatCreateTaoMF(),MatTaoMFGetH(),
          MatTaoMFSetHHistory(), MatTaoMFResetHHistory(),
          MatTaoMFKSPMonitor(), TaoSetGradient()
@*/
int MatTaoMFSetGradient(Mat mat,Vec v,int (*func)(TAO_SOLVER,Vec,Vec,void *),void *funcctx)
{
  MatTaoMFCtx ctx;
  int          info;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  if (ctx) {
    ctx->func    = func;
    ctx->funcctx = funcctx;
    ctx->funcvec = v;
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFSetGradientError"
/*@
   MatTaoMFSetGradientError - Sets the relative error in gradient
   computations for the approximation of Hessian-vector products using
   finite differences of gradients.

   Collective on Mat

   Input Parameters:
+  mat - the matrix free matrix created via MatCreateTaoMF()
-  error_rel - relative error (should be set to the square root of
               the relative error in the gradient evaluations)

   Options Database Keys:
+  -tao_mf_err <error_rel> - Sets error_rel

   Level: advanced

   Notes:
   The default matrix-free matrix-vector product routine computes
.vb
     H(u)*a = [g(u+h*a) - g(u)] / h,
     where
     g(u) denotes the gradient of f(u),
     H(u) denotes the Hessian of f(u), and

     h = error_rel*u'a/||a||^2                        if  |u'a| > umin*||a||_{1}
       = error_rel*umin*sign(u'a)*||a||_{1}/||a||^2   otherwise
 where
     error_rel = square root of relative error in gradient evaluation
     umin = minimum iterate parameter
.ve

.keywords: TAO_SOLVER, matrix-free, parameters

.seealso: MatCreateTaoMF(),MatTaoMFGetH(),
          MatTaoMFSetHHistory(), MatTaoMFResetHHistory(),
          MatTaoMFKSPMonitor()
@*/
int MatTaoMFSetGradientError(Mat mat,double error)
{
  MatTaoMFCtx ctx;
  int          info;

  PetscFunctionBegin;
  info = MatShellGetContext(mat,(void **)&ctx);CHKERRQ(info);
  if (ctx) {
    if (error != TAO_DEFAULT) ctx->error_rel = error;
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFSetHHistory"
/*@
   MatTaoMFSetHHistory - Sets an array to collect a history of the
   differencing values (h) computed for the matrix-free product.

   Collective on Mat 

   Input Parameters:
+  J - the matrix-free matrix context
.  histroy - space to hold the history
-  nhistory - number of entries in history, if more entries are generated than
              nhistory, then the later ones are discarded

   Level: advanced

   Notes:
   Use MatTaoMFResetHHistory() to reset the history counter and collect
   a new batch of differencing parameters, h.

.keywords: TAO_SOLVER, matrix-free, h history, differencing history

.seealso: MatTaoMFGetH(), MatCreateTaoMF(),
          MatTaoMFResetHHistory(),
          MatTaoMFKSPMonitor(), MatTaoMFSetGradientError()

@*/
int MatTaoMFSetHHistory(Mat J,Scalar *history,int nhistory)
{
  int          info;
  MatTaoMFCtx ctx;

  PetscFunctionBegin;

  info = MatShellGetContext(J,(void **)&ctx);CHKERRQ(info);
  /* no context indicates that it is not the "matrix free" matrix type */
  if (!ctx) PetscFunctionReturn(0);
  ctx->historyh    = history;
  ctx->maxcurrenth = nhistory;
  ctx->currenth    = 0;

  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "MatTaoMFResetHHistory"
/*@
   MatTaoMFResetHHistory - Resets the counter to zero to begin 
   collecting a new set of differencing histories.

   Collective on Mat 

   Input Parameters:
.  J - the matrix-free matrix context

   Level: advanced

   Notes:
   Use MatTaoMFSetHHistory() to create the original history counter.

.keywords: TAO_SOLVER, matrix-free, h history, differencing history

.seealso: MatTaoMFGetH(), MatCreateTaoMF(),
          MatTaoMFSetHHistory(),
          MatTaoMFKSPMonitor(), MatTaoMFSetGradientError()

@*/
int MatTaoMFResetHHistory(Mat J)
{
  int          info;
  MatTaoMFCtx ctx;

  PetscFunctionBegin;

  info = MatShellGetContext(J,(void **)&ctx);CHKERRQ(info);
  /* no context indicates that it is not the "matrix free" matrix type */
  if (!ctx) PetscFunctionReturn(0);
  ctx->ncurrenth    = 0;

  PetscFunctionReturn(0);
}

