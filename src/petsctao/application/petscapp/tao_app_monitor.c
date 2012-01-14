#include "tao_app_impl.h"     /*I  "tao.h"  I*/
extern int TAO_APP_COOKIE;

#undef __FUNCT__  
#define __FUNCT__ "TaoAppCheckConvergence"
/*@C
   TaoAppCheckConvergence - Check convergence of application.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  GL - the gradient of the Lagrangian function

   Output Parameters:
.  flag - set to PETSC_TRUE if the TAO solver should stop.

   Level: developer

.keywords: TAO_APPLICATION, monitor

@*/
int TaoAppCheckConvergence(TAO_APPLICATION taoapp, Vec GL, PetscTruth *flag){
  int info;
  PetscInt n;
  double glnorm;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  info=VecNorm(GL,NORM_2,&glnorm);  CHKERRQ(info);
  info=VecGetSize(GL,&n);  CHKERRQ(info);
  *flag=PETSC_FALSE;
  if (glnorm/sqrt((double)n) < taoapp->grtol) *flag=PETSC_TRUE;
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetRelativeTolerance"
/*@C
   TaoAppSetRelativeTolerance - Set convergence tolerance

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
.  grtol - relative tolerance

   Level: intermediate

.keywords: TAO_APPLICATION, monitor

@*/
int TaoAppSetRelativeTolerance(TAO_APPLICATION taoapp, double grtol){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  taoapp->grtol=grtol;
  PetscFunctionReturn(0);
}

