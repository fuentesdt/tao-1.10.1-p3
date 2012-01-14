#include "tao_app_impl.h"       /*I  "tao.h"  I*/

extern int Tao_JacobianEval, Tao_FunctionEval;
extern int TAO_APP_COOKIE;


#undef __FUNCT__  
#define __FUNCT__ "TaoAppComputeVariableBounds"
/*@C
   TaoAppComputeGradient - Compute the gradient of the objective function using the
   routine set by TaoApplicationSetGradientRoutine().

   Collective on TAO_APPLICATION

   Input Parameters:
+  taopp - the TAO_APPLICATION context
-  X - the point where the objective should be evaluated

   Output Parameter:
.  f - function value

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeGradient(), TaoAppComputeObjectiveAndGradient()
@*/
int TaoAppComputeVariableBounds(TAO_APPLICATION taoapp, Vec XL, Vec XU){
  int info;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(XL,VEC_COOKIE,2);
  PetscValidHeaderSpecific(XU,VEC_COOKIE,3);
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  info = PetscInfo(XL,"TAO:  Compute variable bounds of TAO_APPLICATION object.\n"); CHKERRQ(info);
  if (taoapp->computebounds){
    info = (*taoapp->computebounds)(taoapp,XL,XU,taoapp->boundctx); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetVariableBoundsRoutine"
/*@C
   TaoAppSetVariableBoundsRoutine - Sets a routine that evaluates the function at
the specified point.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  func - bound evaluation routine
-  ctx - [optional] user-defined context for private data for the 
         bounds evaluation routine (may be TAO_NULL)

   Calling sequence of func:
$     func (Vec xl,Vec xu, void *ctx);

+  tao - the TAO_APPLICATION context
.  xl - lower bound vector
.  xu - upper bound vector
-  ctx - [optional] user-defined function context 

   Level: beginner

.seealso: TaoGetVariableBoundVecs()

.keywords: TAO_APPLICATION, set, bounds
@*/
int TaoAppSetVariableBoundsRoutine(TAO_APPLICATION taoapp, int (*func)(TAO_APPLICATION,Vec,Vec,void*),void *ctx){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->computebounds=func;
  taoapp->boundctx=ctx;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetFunctionVec"
/*@
   TaoAppSetFunctionVec - Set the Vec that will be used to store the constraint function.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
-  r - vector to constrainf function values

   Level: intermediate

.keywords: TAO_APPLICATION, set, function

.seealso: TaoAppSetJacobianMat(), TaoAppSetConstraintRoutine()

@*/
int TaoAppSetFunctionVec(TAO_APPLICATION taoapp, Vec r){
  int info;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (r){
    PetscValidHeaderSpecific(r,VEC_COOKIE,2);
    PetscObjectReference((PetscObject)r);
  }
  if (taoapp->R){
    info=VecDestroy(taoapp->R); CHKERRQ(info);
  }
  taoapp->R=r;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetFunctionVec"
/*@C
   TaoAppGetFunctionVec - Get the Vec that used to store the constraint function.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
-  r - vector to constrainf function values

   Level: intermediate

.keywords: TAO_APPLICATION, set, function

.seealso: TaoAppSetJacobianMat(), TaoAppSetConstraintRoutine()

@*/
int TaoAppGetFunctionVec(TAO_APPLICATION taoapp, Vec *r){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (r){
    *r=taoapp->R;
  }
  PetscFunctionReturn(0);
}



#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetConstraintRoutine"
/*@C
   TaoAppSetConstraintRoutine - Sets the routine that evaluates
   equality constraint functions.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  func - the constraint function evaluation routine
-  ctx - [optional] user-defined function context 

   Calling sequence of func:
$     func (TAO_APPLICATION taoapp,Vec x,Vec r,void *ctx);

+  taoapp - the TAO_APPLICATION  context
.  x - input vector
.  r - constraint vector
-  ctx - user-defined function gradient context set from TaoAppSetConstraintRoutine()

   Level: intermediate

.keywords: TAO_APPLICATION, set, function

.seealso: TaoAppSetFunctionVec(), TaoAppSetJacobianRoutine()

@*/
int TaoAppSetConstraintRoutine(TAO_APPLICATION taoapp, int (*func)(TAO_APPLICATION,Vec,Vec,void*),void *ctx){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->usrfvctx = ctx;
  taoapp->computevfunc = func;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppComputeFunction"
/*@C
   TaoAppComputeFunction - Compute the constraint vector using the
   routine set by TaoApplicationSetConstraintsRoutine().

   Collective on TAO_APPLICATION

   Input Parameters:
+  taopp - the TAO_APPLICATION context
-  X - the point where the objective should be evaluated

   Output Parameter:
.  R - constraint vector

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeJacobian() TaoAppSetConstraintRoutine()
@*/
int TaoAppComputeFunction(TAO_APPLICATION taoapp, Vec X, Vec R){
  int     info;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  PetscValidHeaderSpecific(X,VEC_COOKIE,2);
  PetscValidHeaderSpecific(R,VEC_COOKIE,3);
  if (taoapp->computevfunc){
    info = PetscLogEventBegin(Tao_FunctionEval,taoapp,X,R,0);CHKERRQ(info);
    PetscStackPush("Tao user ConstraintsFunction");
    info = (*taoapp->computevfunc)(taoapp,X,R,taoapp->usrfvctx);
    CHKERRQ(info);
    PetscStackPop;
    info = PetscLogEventEnd(Tao_FunctionEval,taoapp,X,R,0);
  } else {
    SETERRQ(1,"TAO ERROR: Must set Constraint function");
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetJacobianMat"
/*@
   TaoAppSetJacobianMat - Sets the matrix to be used for the Jacobian.

   Collective on TAO_APPLICATION and Mat

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  J - Jacobian matrix
-  JP - Preconditioner for Jacobian

   Level: intermediate

.keywords: TAO_APPLICATION, Jacobian

.seealso: TaoAppSetJacobianRoutine(), TaoAppSetFunctionVec()
@*/
int TaoAppSetJacobianMat(TAO_APPLICATION taoapp,Mat J, Mat JP){
  int info;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (J || JP){
    PetscValidHeaderSpecific(J,MAT_COOKIE,2);
    PetscValidHeaderSpecific(JP,MAT_COOKIE,3);
    PetscObjectReference((PetscObject)J);
    PetscObjectReference((PetscObject)JP);
  }

  if (taoapp->J || taoapp->JP){
    info=MatDestroy(taoapp->J);CHKERRQ(info);
    info=MatDestroy(taoapp->JP);CHKERRQ(info);
  }
  taoapp->J=J;
  taoapp->JP=JP;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetJacobianMat"
/*@
   TaoAppGetJacobianMat - Get the matrix to be used for the Jacobian.

   Collective on TAO_APPLICATION and Mat

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
-  J - Jacobian matrix

   Level: intermediate

.keywords: TAO_APPLICATION, Jacobian

.seealso: TaoAppSetJacobianMat()
@*/
int TaoAppGetJacobianMat(TAO_APPLICATION taoapp,Mat *J, Mat *JP){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (J) { *J=taoapp->J;}
  if (JP){ *JP=taoapp->JP; }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetJacobianRoutine"
/*@C
   TaoAppSetJacobianRoutine - Sets the function to compute the Jacobian of
the equality constraint function as well as the
   location to store the matrix.

   Collective on TAO_APPLICATION and Mat

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  jac - Jacobian evaluation routine
-  ctx - [optional] user-defined context for private data for the 
         Hessian evaluation routine (may be TAO_NULL)

   Calling sequence of func:
$    jac (TAO_APPLICATION taoapp,Vec x,Mat *J,Mat *JPre, void *ctx);

+  taoapp - the TAO_APPLICATION  context
.  x - input vector
.  J - Jacobian matrix
.  JPre - matrix to precondition J
-  ctx - [optional] user-defined Hessian context

   Notes: 

   The function jac() takes Mat * as the matrix arguments rather than Mat.  
   This allows the Jacobian evaluation routine to replace J with a 
   completely new new matrix structure (not just different matrix elements)
   when appropriate, for instance, if the nonzero structure is changing
   throughout the global iterations.

   Level: intermediate

.keywords: TAO_APPLICATION, Jacobian

.seealso: TaoAppSetJacobianMat(), TaoAppSetConstraintRoutine()
@*/
int TaoAppSetJacobianRoutine(TAO_APPLICATION taoapp,int (*jac)(TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure*,void*),void *ctx){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->computejacobian=jac;
  taoapp->usrjctx=ctx;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppComputeJacobian"
/*@C
   TaoAppComputeJacobian - Compute the Jacobian of the nonlinear equations using the
   routine set by TaoApplicationSetJacobianRoutine().

   Collective on TAO_APPLICATION

   Input Parameters:
+  taopp - the TAO_APPLICATION context
.  X - the variable vector
.  H - the Jacobian matrix
.  HP - the preconditioner for the Jacobian matrix.
-  flag - flag used in KSPSetOperators()

   Output Parameter:
+  H - the Jacobian matrix
.  HP - the preconditioner for the Jacobian matrix.
-  flag - flag used in KSPSetOperators()

   Level: developer

.keywords: TAO_APPLICATION, objective

.seealso: TaoAppComputeFunction(), TaoAppSetJacobianRoutine()
@*/
int TaoAppComputeJacobian(TAO_APPLICATION taoapp, Vec X, Mat *JJ, Mat *JJPre, MatStructure*flag){

  int     info;
  Mat J=*JJ,JPre=*JJPre;
  MatStructure pflag=*flag;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,VEC_COOKIE,2);
  PetscValidHeaderSpecific(J,MAT_COOKIE,3);
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (taoapp->computejacobian){
    PetscStackPush("TAO User Jacobian Evaluation");
    info = PetscLogEventBegin(Tao_JacobianEval,taoapp,X,J,0);CHKERRQ(info);
    info = (*taoapp->computejacobian)(taoapp,X,&J,&JPre, &pflag, taoapp->usrjctx);
    CHKERRQ(info);
    info = PetscLogEventEnd(Tao_JacobianEval,taoapp,X,J,0);CHKERRQ(info);
    PetscStackPop;
  } else {
    SETERRQ(1,"TAO Error:  No Jacobian Routine Available.");
  }
  *JJ=J;*JJPre=JPre; *flag=pflag;
  PetscFunctionReturn(0);
}

