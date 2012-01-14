/*$Id$*/

#include "src/tao_impl.h"    /*I  "tao.h"  I*/

#undef __FUNC__  
#define __FUNC__ "TaoDefaultComputeHessian"
/*@C
   TaoDefaultComputeHessian - Computes the Hessian using finite differences. 

   Collective on TAO_SOLVER

   Input Parameters:
+  x1 - compute Hessian at this point
-  ctx - application's gradient context, as set with TaoSetGradient()

   Output Parameters:
+  H - Hessian matrix (not altered in this routine)
.  B - newly computed Hessian matrix to use with preconditioner (generally the same as H)
-  flag - flag indicating whether the matrix sparsity structure has changed

   Options Database Key:
$  -tao_fd - Activates TaoDefaultComputeHessian()

   Level: advanced

   Notes:
   This routine is slow and expensive, and is not currently optimized
   to take advantage of sparsity in the problem.  Although
   TaoDefaultComputeHessian() is not recommended for general use
   in large-scale applications, It can be useful in checking the
   correctness of a user-provided Hessian.

Concepts: TAO_SOLVER, finite differences, Hessian

.seealso: TaoSetHessian(), TaoDefaultComputeHessianColor()
@*/
int TaoDefaultComputeHessian(TAO_SOLVER tao,Vec x1,Mat *H,Mat *B,
                              MatStructure *flag,void *ctx)
{
  Vec      j1a,j2a,x2;
  int      i,info,N,start,end,j;
  Scalar   dx, mone = -1.0,*y,scale,*xx,wscale;
  double   amax, epsilon = 1.e-8; /* assumes double precision */
  double   dx_min = 1.e-16, dx_par = 1.e-1;
  MPI_Comm comm;
  int      (*eval_fct)(TAO_SOLVER,Vec,Vec)=0;

  PetscFunctionBegin;
  eval_fct = TaoComputeGradient;

  info = PetscObjectGetComm((PetscObject)x1,&comm);CHKERRQ(info);
  info = MatZeroEntries(*B);CHKERRQ(info);
  if (!tao->nvwork) {
    info = VecDuplicateVecs(x1,3,&tao->vwork);CHKERRQ(info);
    tao->nvwork = 3;
    PLogObjectParents(tao,3,tao->vwork);
  }
  j1a = tao->vwork[0]; j2a = tao->vwork[1]; x2 = tao->vwork[2];

  info = VecGetSize(x1,&N);CHKERRQ(info);
  info = VecGetOwnershipRange(x1,&start,&end);CHKERRQ(info);
  info = eval_fct(tao,x1,j1a);CHKERRQ(info);

  /* Compute Hessian approximation, 1 column at a time. 
      x1 = current iterate, j1a = F(x1)
      x2 = perturbed iterate, j2a = F(x2)
   */
  for ( i=0; i<N; i++ ) {
    info = VecCopy(x1,x2);CHKERRQ(info);
    if ( i>= start && i<end) {
      info = VecGetArray(x1,&xx);CHKERRQ(info);
      dx = xx[i-start];
      info = VecRestoreArray(x1,&xx);CHKERRQ(info);
#if !defined(PETSC_USE_COMPLEX)
      if (dx < dx_min && dx >= 0.0) dx = dx_par;
      else if (dx < 0.0 && dx > -dx_min) dx = -dx_par;
#else
      if (PetscAbsScalar(dx) < dx_min && PetscReal(dx) >= 0.0) dx = dx_par;
      else if (PetscReal(dx) < 0.0 && PetscAbsScalar(dx) < dx_min) dx = -dx_par;
#endif
      dx *= epsilon;
      wscale = 1.0/dx;
      info = VecSetValues(x2,1,&i,&dx,ADD_VALUES);CHKERRQ(info);
    } else {
      wscale = 0.0;
    }
    info = eval_fct(tao,x2,j2a);CHKERRQ(info);
    info = VecAXPY(j2a, mone, j1a);CHKERRQ(info);
    /* Communicate scale to all processors */
#if !defined(PETSC_USE_COMPLEX)
    info = MPI_Allreduce(&wscale,&scale,1,MPI_DOUBLE,MPI_SUM,comm);CHKERRQ(info);
#else
    info = MPI_Allreduce(&wscale,&scale,2,MPI_DOUBLE,MPI_SUM,comm);CHKERRQ(info);
#endif
    info = VecScale(j2a, scale);CHKERRQ(info);
    info = VecGetArray(j2a,&y);CHKERRQ(info);
    info = VecNorm(j2a,NORM_INFINITY,&amax);CHKERRQ(info); amax *= 1.e-14;
    for ( j=start; j<end; j++ ) {
      if (PetscAbsScalar(y[j-start]) > amax) {
        info = MatSetValues(*B,1,&j,1,&i,y+j-start,INSERT_VALUES);CHKERRQ(info);
      }
    }
    info = VecRestoreArray(j2a,&y);CHKERRQ(info);
  }
  info  = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(info);
  info  = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(info);
  *flag =  DIFFERENT_NONZERO_PATTERN;
  PetscFunctionReturn(0);
}
