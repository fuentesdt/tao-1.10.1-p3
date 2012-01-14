#include "tao_general.h"
#include "private/vecimpl.h"    /*I "tao_solver.h" I*/
#include "petscksp.h"



#undef __FUNCT__  
#define __FUNCT__ "VecCompare"
int  VecCompare(Vec V1,Vec V2, PetscTruth *flg){
  int info;
  PetscInt n1,n2,N1,N2;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V1,VEC_COOKIE,1); 
  PetscValidHeaderSpecific(V2,VEC_COOKIE,2); 
  info = VecGetSize(V1,&N1);CHKERRQ(info);
  info = VecGetSize(V2,&N2);CHKERRQ(info);
  info = VecGetLocalSize(V1,&n1);CHKERRQ(info);
  info = VecGetLocalSize(V2,&n2);CHKERRQ(info);
  if (N1==N2 && n1==n2) 
    *flg=PETSC_TRUE;
  else
    *flg=PETSC_FALSE;

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "VecMedian"
int VecMedian(Vec Vec1, Vec Vec2, Vec Vec3, Vec VMedian)
{
  int i,info;
  PetscInt n,low1,low2,low3,low4,high1,high2,high3,high4;
  PetscScalar *v1,*v2,*v3,*vmed;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Vec1,VEC_COOKIE,1); 
  PetscValidHeaderSpecific(Vec2,VEC_COOKIE,2); 
  PetscValidHeaderSpecific(Vec3,VEC_COOKIE,3); 
  PetscValidHeaderSpecific(VMedian,VEC_COOKIE,4); 

  if (Vec1==Vec2 || Vec1==Vec3){
    info=VecCopy(Vec1,VMedian); CHKERRQ(info); 
    PetscFunctionReturn(0);
  }
  if (Vec2==Vec3){
    info=VecCopy(Vec2,VMedian); CHKERRQ(info); 
    PetscFunctionReturn(0);
  }

  PetscValidType(Vec1,1);
  PetscValidType(Vec2,2);
  PetscValidType(VMedian,4);
  PetscCheckSameType(Vec1,1,Vec2,2); PetscCheckSameType(Vec1,1,VMedian,4);
  PetscCheckSameComm(Vec1,1,Vec2,2); PetscCheckSameComm(Vec1,1,VMedian,4);

  info = VecGetOwnershipRange(Vec1, &low1, &high1); CHKERRQ(info);
  info = VecGetOwnershipRange(Vec2, &low2, &high2); CHKERRQ(info);
  info = VecGetOwnershipRange(Vec3, &low3, &high3); CHKERRQ(info);
  info = VecGetOwnershipRange(VMedian, &low4, &high4); CHKERRQ(info);
  if ( low1!= low2 || low1!= low3 || low1!= low4 ||
       high1!= high2 || high1!= high3 || high1!= high4)
    SETERRQ(1,"InCompatible vector local lengths");

  info = VecGetArray(Vec1,&v1); CHKERRQ(info);
  info = VecGetArray(Vec2,&v2); CHKERRQ(info);
  info = VecGetArray(Vec3,&v3); CHKERRQ(info);

  if ( VMedian != Vec1 && VMedian != Vec2 && VMedian != Vec3){
    info = VecGetArray(VMedian,&vmed); CHKERRQ(info);
  } else if ( VMedian==Vec1 ){
    vmed=v1;
  } else if ( VMedian==Vec2 ){
    vmed=v2;
  } else {
    vmed=v3;
  }

  info=VecGetLocalSize(Vec1,&n); CHKERRQ(info);

  for (i=0;i<n;i++){
    vmed[i]=PetscMax(PetscMax(PetscMin(v1[i],v2[i]),PetscMin(v1[i],v3[i])),PetscMin(v2[i],v3[i]));
  }

  info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
  info = VecRestoreArray(Vec2,&v2); CHKERRQ(info);
  info = VecRestoreArray(Vec3,&v2); CHKERRQ(info);
  
  if (VMedian!=Vec1 && VMedian != Vec2 && VMedian != Vec3){
    info = VecRestoreArray(VMedian,&vmed); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------- 
#undef __FUNCT__  
#define __FUNCT__ "VecPointwiseMin"
int VecPointwiseMin(Vec Vec1, Vec Vec2, Vec VMin)
{
  int i,n,low1,low2,low3,high1,high2,high3,info;
  PetscScalar *v1,*v2,*vmin;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Vec1,VEC_COOKIE,1); 
  PetscValidHeaderSpecific(Vec2,VEC_COOKIE,2); 
  PetscValidHeaderSpecific(VMin,VEC_COOKIE,3); 

  if (Vec1==Vec2){
    info=VecCopy(Vec1,VMin); CHKERRQ(info); 
    PetscFunctionReturn(0);
  }

  PetscCheckSameType(Vec1,1,Vec2,2); PetscCheckSameType(Vec1,1,VMin,3);
  PetscCheckSameComm(Vec1,1,Vec2,2); PetscCheckSameComm(Vec1,1,VMin,3);

  info = VecGetOwnershipRange(Vec1, &low1, &high1); CHKERRQ(info);
  info = VecGetOwnershipRange(Vec2, &low2, &high2); CHKERRQ(info);
  info = VecGetOwnershipRange(VMin, &low3, &high3); CHKERRQ(info);

  if ( low1!= low2 || low1!= low3 || high1!= high2 || high1!= high3)
    SETERRQ(1,"InCompatible vector local lengths");

  info = VecGetArray(Vec1,&v1); CHKERRQ(info);
  info = VecGetArray(Vec2,&v2); CHKERRQ(info);

  if ( VMin != Vec1 && VMin != Vec2){
    info = VecGetArray(VMin,&vmin); CHKERRQ(info);
  } else if ( VMin==Vec1 ){
    vmin=v1;
  } else {
    vmin =v2;
  }

  info=VecGetLocalSize(Vec1,&n); CHKERRQ(info);

  for (i=0; i<n; i++){
    vmin[i] = PetscMin(v1[i],v2[i]);
  }

  info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
  info = VecRestoreArray(Vec2,&v2); CHKERRQ(info);

  if (VMin!=Vec1 && VMin != Vec2){
    info = VecRestoreArray(VMin,&vmin); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}
*/
/* ---------------------------------------------------------- 
#undef __FUNCT__  
#define __FUNCT__ "VecPointwiseMax"
int VecPointwiseMax(Vec Vec1, Vec Vec2, Vec VMax)
{
  int i,n,low1,low2,low3,high1,high2,high3,info;
  PetscScalar *v1,*v2,*vmax;

  PetscFunctionBegin;

  if (Vec1==Vec2){
    info=VecCopy(Vec1,VMax); CHKERRQ(info); 
    PetscFunctionReturn(0);
  }

  PetscValidHeaderSpecific(Vec1,VEC_COOKIE,1); 
  PetscValidHeaderSpecific(Vec2,VEC_COOKIE,2); 
  PetscValidHeaderSpecific(VMax,VEC_COOKIE,3); 
  PetscValidType(VMax,4);
  PetscCheckSameType(Vec1,1,Vec2,2); PetscCheckSameType(Vec1,1,VMax,3);
  PetscCheckSameComm(Vec1,1,Vec2,2); PetscCheckSameComm(Vec1,1,VMax,3);

  info = VecGetOwnershipRange(Vec1, &low1, &high1); CHKERRQ(info);
  info = VecGetOwnershipRange(Vec2, &low2, &high2); CHKERRQ(info);
  info = VecGetOwnershipRange(VMax, &low3, &high3); CHKERRQ(info);

  if ( low1!= low2 || low1!= low3 || high1!= high2 || high1!= high3)
    SETERRQ(1,"Vectors must be identically loaded over processors");

  info = VecGetArray(Vec1,&v1); CHKERRQ(info);
  info = VecGetArray(Vec2,&v2); CHKERRQ(info);

  if ( VMax != Vec1 && VMax != Vec2){
    info = VecGetArray(VMax,&vmax); CHKERRQ(info);
  } else if ( VMax==Vec1 ){ vmax=v1;
  } else { vmax =v2;
  }

  info=VecGetLocalSize(Vec1,&n); CHKERRQ(info);

  for (i=0; i<n; i++){
    vmax[i] = PetscMax(v1[i],v2[i]);
  }

  info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
  info = VecRestoreArray(Vec2,&v2); CHKERRQ(info);

  if ( VMax != Vec1 && VMax != Vec2){
    info = VecRestoreArray(VMax,&vmax); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}*/

#undef __FUNCT__
#define __FUNCT__ "Fischer"
inline static PetscScalar Fischer(PetscScalar a, PetscScalar b)
{
   // Method suggested by Bob Vanderbei
   if (a + b <= 0) {
     return sqrt(a*a + b*b) - (a + b);
   }
   return -2.0*a*b / (sqrt(a*a + b*b) + (a + b));
}

#undef __FUNCT__  
#define __FUNCT__ "VecFischer"
int VecFischer(Vec X, Vec F, Vec L, Vec U, Vec FF)
{
  PetscScalar *x, *f, *l, *u, *ff;
  PetscScalar xval, fval, lval, uval;
  int info, i;
  PetscInt low[5], high[5], n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X, VEC_COOKIE,1); 
  PetscValidHeaderSpecific(F, VEC_COOKIE,2); 
  PetscValidHeaderSpecific(L, VEC_COOKIE,3); 
  PetscValidHeaderSpecific(U, VEC_COOKIE,4); 
  PetscValidHeaderSpecific(FF, VEC_COOKIE,4); 

  info = VecGetOwnershipRange(X, low, high); CHKERRQ(info);
  info = VecGetOwnershipRange(F, low + 1, high + 1); CHKERRQ(info);
  info = VecGetOwnershipRange(L, low + 2, high + 2); CHKERRQ(info);
  info = VecGetOwnershipRange(U, low + 3, high + 3); CHKERRQ(info);
  info = VecGetOwnershipRange(FF, low + 4, high + 4); CHKERRQ(info);

  for (i = 1; i < 4; ++i) {
    if (low[0] != low[i] || high[0] != high[i])
      SETERRQ(1,"Vectors must be identically loaded over processors");
  }

  info = VecGetArray(X, &x); CHKERRQ(info);
  info = VecGetArray(F, &f); CHKERRQ(info);
  info = VecGetArray(L, &l); CHKERRQ(info);
  info = VecGetArray(U, &u); CHKERRQ(info);
  info = VecGetArray(FF, &ff); CHKERRQ(info);

  info = VecGetLocalSize(X, &n); CHKERRQ(info);

  for (i = 0; i < n; ++i) {
    xval = x[i]; fval = f[i];
    lval = l[i]; uval = u[i];

    if ((lval <= -TAO_INFINITY) && (uval >= TAO_INFINITY)) {
      ff[i] = -fval;
    } 
    else if (lval <= -TAO_INFINITY) {
      ff[i] = -Fischer(uval - xval, -fval);
    } 
    else if (uval >=  TAO_INFINITY) {
      ff[i] =  Fischer(xval - lval,  fval);
    } 
    else if (lval == uval) {
      ff[i] = lval - xval;
    }
    else {
      fval  =  Fischer(uval - xval, -fval);
      ff[i] =  Fischer(xval - lval,  fval);
    }
  }
  
  info = VecRestoreArray(X, &x); CHKERRQ(info);
  info = VecRestoreArray(F, &f); CHKERRQ(info);
  info = VecRestoreArray(L, &l); CHKERRQ(info);
  info = VecRestoreArray(U, &u); CHKERRQ(info);
  info = VecRestoreArray(FF, &ff); CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SFischer"
inline static PetscScalar SFischer(PetscScalar a, PetscScalar b, PetscScalar c)
{
   // Method suggested by Bob Vanderbei
   if (a + b <= 0) {
     return sqrt(a*a + b*b + 2.0*c*c) - (a + b);
   }
   return 2.0*(c*c - a*b) / (sqrt(a*a + b*b + 2.0*c*c) + (a + b));
}

#undef __FUNCT__
#define __FUNCT__ "VecSFischer"
int VecSFischer(Vec X, Vec F, Vec L, Vec U, double mu, Vec FF)
{
  PetscScalar *x, *f, *l, *u, *ff;
  PetscScalar xval, fval, lval, uval;
  int info, i;
  PetscInt low[5], high[5], n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X, VEC_COOKIE,1);
  PetscValidHeaderSpecific(F, VEC_COOKIE,2);
  PetscValidHeaderSpecific(L, VEC_COOKIE,3);
  PetscValidHeaderSpecific(U, VEC_COOKIE,4);
  PetscValidHeaderSpecific(FF, VEC_COOKIE,6);

  info = VecGetOwnershipRange(X, low, high); CHKERRQ(info);
  info = VecGetOwnershipRange(F, low + 1, high + 1); CHKERRQ(info);
  info = VecGetOwnershipRange(L, low + 2, high + 2); CHKERRQ(info);
  info = VecGetOwnershipRange(U, low + 3, high + 3); CHKERRQ(info);
  info = VecGetOwnershipRange(FF, low + 4, high + 4); CHKERRQ(info);

  for (i = 1; i < 4; ++i) {
    if (low[0] != low[i] || high[0] != high[i])
      SETERRQ(1,"Vectors must be identically loaded over processors");
  }

  info = VecGetArray(X, &x); CHKERRQ(info);
  info = VecGetArray(F, &f); CHKERRQ(info);
  info = VecGetArray(L, &l); CHKERRQ(info);
  info = VecGetArray(U, &u); CHKERRQ(info);
  info = VecGetArray(FF, &ff); CHKERRQ(info);

  info = VecGetLocalSize(X, &n); CHKERRQ(info);

  for (i = 0; i < n; ++i) {
    xval = (*x++); fval = (*f++);
    lval = (*l++); uval = (*u++);

    if ((lval <= -TAO_INFINITY) && (uval >= TAO_INFINITY)) {
      (*ff++) = -fval - mu*xval;
    } 
    else if (lval <= -TAO_INFINITY) {
      (*ff++) = -SFischer(uval - xval, -fval, mu);
    } 
    else if (uval >=  TAO_INFINITY) {
      (*ff++) =  SFischer(xval - lval,  fval, mu);
    } 
    else if (lval == uval) {
      (*ff++) = lval - xval;
    } 
    else {
      fval    =  SFischer(uval - xval, -fval, mu);
      (*ff++) =  SFischer(xval - lval,  fval, mu);
    }
  }
  x -= n; f -= n; l -=n; u -= n; ff -= n;

  info = VecRestoreArray(X, &x); CHKERRQ(info);
  info = VecRestoreArray(F, &f); CHKERRQ(info);
  info = VecRestoreArray(L, &l); CHKERRQ(info);
  info = VecRestoreArray(U, &u); CHKERRQ(info);
  info = VecRestoreArray(FF, &ff); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecBoundProjection"
int VecBoundProjection(Vec G, Vec X, Vec XL, Vec XU, Vec GP){

  int i,info;
  PetscInt n;
  PetscScalar *xptr,*xlptr,*xuptr,*gptr,*gpptr;
  PetscScalar xval,gpval;

  /* Project variables at the lower and upper bound */

  PetscFunctionBegin;
  info = VecGetLocalSize(X,&n); CHKERRQ(info);

  info=VecGetArray(X,&xptr); CHKERRQ(info);
  info=VecGetArray(XL,&xlptr); CHKERRQ(info);
  info=VecGetArray(XU,&xuptr); CHKERRQ(info);
  info=VecGetArray(G,&gptr); CHKERRQ(info);
  if (G!=GP){
    info=VecGetArray(GP,&gpptr); CHKERRQ(info);
  } else { gpptr=gptr; }

  for (i=0; i<n; ++i){
    gpval = gptr[i]; xval = xptr[i]; 

    if (gpval>0 && xval<=xlptr[i]){
      gpval = 0;
    } else if (gpval<0 && xval>=xuptr[i]){
      gpval = 0;
    }
    gpptr[i] = gpval;

    /*
    if (xlptr[i] >= xuptr[i]){
      gpptr[i]=0.0;
    } else if (xptr[i] <= xlptr[i]+eps){
      gpptr[i] = PetscMin(gptr[i],zero);
    } else if (xptr[i] >= xuptr[i]-eps){
      gpptr[i] = PetscMax(gptr[i],zero);
    } else {
      gpptr[i] = gptr[i];
    }
    */
  }

  info=VecRestoreArray(X,&xptr); CHKERRQ(info);
  info=VecRestoreArray(XL,&xlptr); CHKERRQ(info);
  info=VecRestoreArray(XU,&xuptr); CHKERRQ(info);
  info=VecRestoreArray(G,&gptr); CHKERRQ(info);
  if (G!=GP){
    info=VecRestoreArray(GP,&gpptr); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecISAXPY"
int VecISAXPY(Vec Y, PetscScalar alpha, Vec X, IS YIS){

  int i,info;
  PetscInt in1,in2,xlow,xhigh,ylow,yhigh;
  const PetscInt *yis;
  PetscScalar *x,*y;

  PetscFunctionBegin;
  info=ISGetLocalSize(YIS,&in1); CHKERRQ(info);
  info=VecGetLocalSize(X,&in2); CHKERRQ(info);

  if ( in1 != in2 )
    SETERRQ(1,"Index set and X vector must be identically loaded over processors");

  info=VecGetOwnershipRange(X, &xlow, &xhigh); CHKERRQ(info);
  info=VecGetOwnershipRange(Y, &ylow, &yhigh); CHKERRQ(info);

  info=ISGetIndices(YIS, &yis); CHKERRQ(info);

  info=VecGetArray(X,&x); CHKERRQ(info);
  if (X != Y){
    info=VecGetArray(Y,&y); CHKERRQ(info);
  } else {
    y=x;
  }

  for (i=0; i<in1; i++){
    if (yis[i]>=ylow && yis[i]<yhigh && i+xlow < xhigh){
      y[yis[i]-ylow] = y[yis[i]-ylow] + (alpha)*x[i];
    } else {
      SETERRQ(1,"IS index out of range");
    }
  }
  
  info=ISRestoreIndices(YIS, &yis); CHKERRQ(info);

  info=VecRestoreArray(X,&x); CHKERRQ(info);
  if (X != Y){
    info=VecRestoreArray(Y,&y); CHKERRQ(info);
  }

  info = PetscLogFlops(2*in1); CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "VecCreateSubVec"
int VecCreateSubVec(Vec A,IS S,Vec *B)
{
  int     i,info;
  PetscInt nloc,n,low,high,nnloc;
  const PetscInt *ptrS;
  PetscScalar  zero=0.0;
  PetscScalar *ptrB,*ptrA;
  const VecType type_name;
  MPI_Comm comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(A,VEC_COOKIE,1);
  PetscValidHeaderSpecific(S,IS_COOKIE,2);
  info=ISGetLocalSize(S,&nloc);CHKERRQ(info);
  info=ISGetSize(S,&n);CHKERRQ(info);
  info=PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(info);
  info=VecCreate(comm,B);CHKERRQ(info);
  info=VecSetSizes(*B,nloc,n);CHKERRQ(info);
  info=VecGetType(A,&type_name);CHKERRQ(info);
  info=VecSetType(*B,type_name);CHKERRQ(info);
  info=VecSet(*B, zero);CHKERRQ(info);
  
  info=VecGetOwnershipRange(A,&low,&high);CHKERRQ(info);
  info=VecGetLocalSize(A,&nnloc);CHKERRQ(info);
  info=VecGetArray(A,&ptrA);CHKERRQ(info);
  info=VecGetArray(*B,&ptrB);CHKERRQ(info);
  info=ISGetIndices(S,&ptrS);CHKERRQ(info);
  for (i=0;i<nloc;i++){ ptrB[i]=ptrA[ptrS[i]-low]; }
  info=VecRestoreArray(A,&ptrA);CHKERRQ(info);
  info=VecRestoreArray(*B,&ptrB);CHKERRQ(info);
  info=ISRestoreIndices(S,&ptrS);CHKERRQ(info);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecStepMax"
int VecStepMax(Vec X, Vec DX, PetscReal*step){
  int i,info;
  PetscInt nn;
  PetscScalar stepmax=1.0e300;
  PetscScalar *xx,*dx;
  double t1,t2;

  PetscFunctionBegin;
  info = VecGetLocalSize(X,&nn);CHKERRQ(info);
  info = VecGetArray(X,&xx);CHKERRQ(info);
  info = VecGetArray(DX,&dx);CHKERRQ(info);
  for (i=0;i<nn;i++){
    if (xx[i] < 0){
      SETERRQ(1,"Vector must be positive");
    } else if (dx[i]<0){ stepmax=PetscMin(stepmax,-xx[i]/dx[i]);
    }
  }
  info = VecRestoreArray(X,&xx);CHKERRQ(info);
  info = VecRestoreArray(DX,&dx);CHKERRQ(info);
  t1=(double)stepmax;
  info = MPI_Allreduce(&t1,&t2,1,MPI_DOUBLE,MPI_MIN,((PetscObject)X)->comm);
  CHKERRQ(info);
  *step=(PetscReal)t2;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatCreateProjectionOperator"
int MatCreateProjectionOperator(Vec A, IS S , Vec B, Mat *MM){

  PetscInt m,n,M,N,low,high;
  int info,i;
  const PetscInt *iptr;
  PetscScalar one=1.0;
  Mat MMM;

  PetscFunctionBegin;  
  PetscValidHeaderSpecific(A,VEC_COOKIE,1);
  PetscValidHeaderSpecific(B,VEC_COOKIE,3);
  PetscValidHeaderSpecific(S,IS_COOKIE,2);
  PetscCheckSameComm(A,1,B,3);

  info = VecGetSize(B,&M);CHKERRQ(info);
  info = VecGetLocalSize(B,&m);CHKERRQ(info);
  info = VecGetSize(A,&N);CHKERRQ(info);
  info = VecGetLocalSize(A,&n);CHKERRQ(info);

  info = MatCreateMPIAIJ(((PetscObject)A)->comm,m,n,M,N,1,PETSC_NULL,1,PETSC_NULL,&MMM);
  CHKERRQ(info);
  /*
  info = MatCreateMPIBDiag(((PetscObject)pv)->comm,m,M,N,M,1,PETSC_NULL,PETSC_NULL,&A);
  */
  
  info = ISGetSize(S,&n);CHKERRQ(info);
  if (n!=M){    
    SETERRQ(1,"Size of index set does not match the second vector.");
  }
  
  info = ISGetLocalSize(S,&n);CHKERRQ(info);
  if (n!=m){    
    SETERRQ(1,"Local size of index set does not match the second vector.");
  }
  info=VecGetOwnershipRange(B,&low,&high);CHKERRQ(info);
  info = ISGetIndices(S,&iptr);CHKERRQ(info);
  for (i=0; i<n; i++){
    info=MatSetValue(MMM,low+i,iptr[i],one,INSERT_VALUES);CHKERRQ(info);
  }
  info = ISRestoreIndices(S,&iptr);CHKERRQ(info);

  info = MatAssemblyBegin(MMM,MAT_FINAL_ASSEMBLY);CHKERRQ(info);
  info = MatAssemblyEnd(MMM,MAT_FINAL_ASSEMBLY);CHKERRQ(info);



  *MM = MMM;

  PetscFunctionReturn(0);  
}

#undef __FUNCT__
#define __FUNCT__ "VecStepMax2"
int VecStepMax2(Vec X, Vec DX, Vec XL, Vec XU, PetscReal *step2){

  int i,info;
  PetscInt nn1,nn2;
  PetscScalar *xx,*dx,*xl,*xu;
  PetscScalar stepmax2=0;
  double t1,t2;

  TaoFunctionBegin;
  info = VecGetArray(X,&xx);CHKERRQ(info);
  info = VecGetArray(XL,&xl);CHKERRQ(info);
  info = VecGetArray(XU,&xu);CHKERRQ(info);
  info = VecGetArray(DX,&dx);CHKERRQ(info);
  info = VecGetLocalSize(X,&nn1);CHKERRQ(info);
  info = VecGetLocalSize(DX,&nn2);CHKERRQ(info);

  for (i=0;i<nn1;i++){
    if (dx[i] > 0){
      stepmax2=PetscMax(stepmax2,(xu[i]-xx[i])/dx[i]);      
    } else if (dx[i]<0){ 
      stepmax2=PetscMax(stepmax2,(xl[i]-xx[i])/dx[i]);
    }
  }
  info = VecRestoreArray(X,&xx);CHKERRQ(info);
  info = VecRestoreArray(XL,&xl);CHKERRQ(info);
  info = VecRestoreArray(XU,&xu);CHKERRQ(info);
  info = VecRestoreArray(DX,&dx);CHKERRQ(info);

  t1=(double)stepmax2;
  info = MPI_Allreduce(&t1,&t2,1,MPI_DOUBLE,MPI_MAX,((PetscObject)X)->comm);
  CHKERRQ(info);
  *step2=(PetscReal)t2;

  TaoFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecBoundStepInfo"
int VecBoundStepInfo(Vec X, Vec XL, Vec XU, Vec DX, PetscReal *boundmin, PetscReal *wolfemin, PetscReal *boundmax){
  int info,i;
  PetscInt n;
  PetscScalar *x,*xl,*xu,*dx;
  double tt,t1=1.0e+20, t2=1.0e+20, t3=0;
  double tt1,tt2,tt3;
  MPI_Comm comm;
  
  PetscFunctionBegin;
  info=VecGetArray(X,&x);CHKERRQ(info);
  info=VecGetArray(XL,&xl);CHKERRQ(info);
  info=VecGetArray(XU,&xu);CHKERRQ(info);
  info=VecGetArray(DX,&dx);CHKERRQ(info);
  info = VecGetLocalSize(X,&n);CHKERRQ(info);
  for (i=0;i<n;i++){
    if (dx[i]>0){
      tt=(xu[i]-x[i])/dx[i];
      t1=PetscMin(t1,tt);
      if (t1>0){
	t2=PetscMin(t2,tt);
      }
      t3=PetscMax(t3,tt);
    } else if (dx[i]<0){
      tt=(xl[i]-x[i])/dx[i];
      t1=PetscMin(t1,tt);
      if (t1>0){
	t2=PetscMin(t2,tt);
      }
      t3=PetscMax(t3,tt);
    }
  }
  info=VecRestoreArray(X,&x);CHKERRQ(info);
  info=VecRestoreArray(XL,&xl);CHKERRQ(info);
  info=VecRestoreArray(XU,&xu);CHKERRQ(info);
  info=VecRestoreArray(DX,&dx);CHKERRQ(info);
  info=PetscObjectGetComm((PetscObject)X,&comm);CHKERRQ(info);
  
  if (boundmin){ info = MPI_Allreduce(&t1,&tt1,1,MPI_DOUBLE,MPI_MIN,comm);CHKERRQ(info);}
  if (wolfemin){ info = MPI_Allreduce(&t2,&tt2,1,MPI_DOUBLE,MPI_MIN,comm);CHKERRQ(info);}
  if (boundmax) { info = MPI_Allreduce(&t3,&tt3,1,MPI_DOUBLE,MPI_MAX,comm);CHKERRQ(info);}

  *boundmin=(PetscReal)tt1;
  *wolfemin=(PetscReal)tt2;
  *boundmax=(PetscReal)tt3;
  info = PetscInfo3(X,"Step Bound Info: Closest Bound: %6.4e, Wolfe: %6.4e, Max: %6.4e \n",*boundmin,*wolfemin,*boundmax); CHKERRQ(info);
  PetscFunctionReturn(0);  
}

#undef __FUNCT__  
#define __FUNCT__ "VecISSetToConstant"
/*@C
   VecISSetToConstant - Sets the elements of a vector, specified by an index set, to a constant

   Input Parameter:
+  S -  a PETSc IS
.  c - the constant
-  V - a Vec

.seealso VecSet()

   Level: advanced
@*/
int VecISSetToConstant(IS S, PetscScalar c, Vec V){
  int info,i;
  PetscInt nloc,low,high;
  const PetscInt *s;
  PetscScalar *v;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,VEC_COOKIE,3); 
  PetscValidHeaderSpecific(S,IS_COOKIE,1); 
  PetscValidType(V,3);
  PetscCheckSameComm(V,3,S,1);

  info = VecGetOwnershipRange(V, &low, &high); CHKERRQ(info);
  info = ISGetLocalSize(S,&nloc);CHKERRQ(info);

  info = ISGetIndices(S, &s); CHKERRQ(info);
  info = VecGetArray(V,&v); CHKERRQ(info);
  for (i=0; i<nloc; i++){
    v[s[i]-low] = c;
    /*
    if (s[i]>=low && s[i]<high){
    v[s[i]-low] = c;
    } else {
      SETERRQ(1,"IS index out of range");
    }
    */
  }
  
  info = ISRestoreIndices(S, &s); CHKERRQ(info);
  info = VecRestoreArray(V,&v); CHKERRQ(info);

  PetscFunctionReturn(0);

}

#undef __FUNCT__  
#define __FUNCT__ "VecPow"
int VecPow(Vec Vec1, double p)
{
  int i,info;
  PetscInt n;
  PetscScalar *v1;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Vec1, VEC_COOKIE, 1); 

  info = VecGetArray(Vec1, &v1); CHKERRQ(info);
  info = VecGetLocalSize(Vec1, &n); CHKERRQ(info);

  if (1.0 == p) {
  }
  else if (-1.0 == p) {
    for (i = 0; i < n; ++i){
      v1[i] = 1.0 / v1[i];
    }
  }
  else if (0.0 == p) {
    for (i = 0; i < n; ++i){
      // Not-a-number left alone
      // Infinity set to one 
      if (v1[i] == v1[i]) {
        v1[i] = 1.0;
      }
    }
  }
  else if (0.5 == p) {
    for (i = 0; i < n; ++i) {
      if (v1[i] >= 0) {
        v1[i] = sqrt(v1[i]);
      }
      else {
        v1[i] = TAO_INFINITY;
      }
    }
  }
  else if (-0.5 == p) {
    for (i = 0; i < n; ++i) {
      if (v1[i] >= 0) {
        v1[i] = 1.0 / sqrt(v1[i]);
      }
      else {
        v1[i] = TAO_INFINITY;
      }
    }
  }
  else if (2.0 == p) {
    for (i = 0; i < n; ++i){
      v1[i] *= v1[i];
    }
  }
  else if (-2.0 == p) {
    for (i = 0; i < n; ++i){
      v1[i] = 1.0 / (v1[i] * v1[i]);
    }
  }
  else {
    for (i = 0; i < n; ++i) {
      if (v1[i] >= 0) {
        v1[i] = pow(v1[i], p);
      }
      else {
        v1[i] = TAO_INFINITY;
      }
    }
  }

  info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
  PetscFunctionReturn(0);
}
