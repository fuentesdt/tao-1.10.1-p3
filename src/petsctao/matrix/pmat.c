#include "tao_general.h"

#ifdef TAO_USE_PETSC
#include "petscmat.h"


int MatD_Fischer(Mat, Vec, Vec, Vec, Vec, Vec, Vec, Vec, Vec);
int MatD_SFischer(Mat, Vec, Vec, Vec, Vec, double, Vec, Vec, Vec, Vec, Vec);

#undef __FUNCT__
#define __FUNCT__ "Fischer"
inline static PetscScalar Fischer(PetscScalar a, PetscScalar b)
{

#ifdef TODD

  if (PetscAbsScalar(a) > PetscAbsScalar(b)) {
    return sqrt(a*a + b*b) - a - b;
  }
  return sqrt(a*a + b*b) - b - a;

#else

   // Method suggested by Bob Vanderbei

   if (a + b <= 0) {
     return sqrt(a*a + b*b) - (a + b);
   }
   return -2.0*a*b / (sqrt(a*a + b*b) + (a + b));

#endif

}

#undef __FUNCT__
#define __FUNCT__ "Norm"
inline static PetscScalar Norm(PetscScalar a, PetscScalar b)
{
  return sqrt(a*a + b*b);
}

#undef __FUNCT__
#define __FUNCT__ "MatD_Fischer"
int MatD_Fischer(Mat m, Vec X, Vec F, Vec L, Vec U,
                 Vec T1, Vec T2, Vec DA, Vec DB)
{
  PetscScalar *x, *f, *l, *u, *da, *db, *t1, *t2;
  PetscScalar ai, bi, ci, di, ei;
  int info, i;
  PetscInt n, low[8], high[8];

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X, VEC_COOKIE,2);
  PetscValidHeaderSpecific(F, VEC_COOKIE,3);
  PetscValidHeaderSpecific(L, VEC_COOKIE,4);
  PetscValidHeaderSpecific(U, VEC_COOKIE,5);
  PetscValidHeaderSpecific(DA, VEC_COOKIE,8);
  PetscValidHeaderSpecific(DB, VEC_COOKIE,9);
  PetscValidHeaderSpecific(T1, VEC_COOKIE,6);
  PetscValidHeaderSpecific(T2, VEC_COOKIE,7);

  info = VecGetOwnershipRange(X, low, high); CHKERRQ(info);
  info = VecGetOwnershipRange(F, low + 1, high + 1); CHKERRQ(info);
  info = VecGetOwnershipRange(L, low + 2, high + 2); CHKERRQ(info);
  info = VecGetOwnershipRange(U, low + 3, high + 3); CHKERRQ(info);
  info = VecGetOwnershipRange(DA, low + 4, high + 4); CHKERRQ(info);
  info = VecGetOwnershipRange(DB, low + 5, high + 5); CHKERRQ(info);
  info = VecGetOwnershipRange(T1, low + 6, high + 6); CHKERRQ(info);
  info = VecGetOwnershipRange(T2, low + 7, high + 7); CHKERRQ(info);

  for (i = 1; i < 8; i++) {
    if (low[0] != low[i] || high[0] != high[i])
      SETERRQ(1,"Vectors must be identically loaded over processors");
  }

  info = VecGetArray(X, &x); CHKERRQ(info);
  info = VecGetArray(F, &f); CHKERRQ(info);
  info = VecGetArray(L, &l); CHKERRQ(info);
  info = VecGetArray(U, &u); CHKERRQ(info);
  info = VecGetArray(DA, &da); CHKERRQ(info);
  info = VecGetArray(DB, &db); CHKERRQ(info);
  info = VecGetArray(T1, &t1); CHKERRQ(info);

  info = VecGetLocalSize(X, &n); CHKERRQ(info);

  for (i = 0; i < n; i++) {
    da[i] = 0;
    db[i] = 0;
    t1[i] = 0;

    if (PetscAbsScalar(f[i]) <= TAO_EPSILON) {
      if (l[i] > -TAO_INFINITY && PetscAbsScalar(x[i] - l[i]) <= TAO_EPSILON) {
        t1[i] = 1;
	da[i] = 1;
      }

      if (u[i] <  TAO_INFINITY && PetscAbsScalar(u[i] - x[i]) <= TAO_EPSILON) {
        t1[i] = 1;
	db[i] = 1;
      }
    }
  }

  info = VecRestoreArray(T1, &t1); CHKERRQ(info);
  info = MatMult(m, T1, T2); CHKERRQ(info);
  info = VecGetArray(T2, &t2); CHKERRQ(info);
  
  for (i = 0; i < n; i++) {
    if ((l[i] <= -TAO_INFINITY) && (u[i] >= TAO_INFINITY)) {
      da[i] = 0;
      db[i] = -1;
    } 
    else if (l[i] <= -TAO_INFINITY) {
      if (db[i] >= 1) {
        ai = Norm(1, t2[i]);

	da[i] = -1/ai - 1;
	db[i] = -t2[i]/ai - 1;
      } 
      else {
        bi = u[i] - x[i];
        ai = Norm(bi, f[i]);
	ai = PetscMax(TAO_EPSILON, ai);

	da[i] = bi / ai - 1;
	db[i] = -f[i] / ai - 1;
      }
    } 
    else if (u[i] >=  TAO_INFINITY) {
      if (da[i] >= 1) {
        ai = Norm(1, t2[i]);

	da[i] = 1 / ai - 1;
	db[i] = t2[i] / ai - 1;
      } 
      else {
        bi = x[i] - l[i];
        ai = Norm(bi, f[i]);
	ai = PetscMax(TAO_EPSILON, ai);

        da[i] = bi / ai - 1;
	db[i] = f[i] / ai - 1;
      }
    } 
    else if (l[i] == u[i]) {
      da[i] = -1;
      db[i] = 0;
    } 
    else {
      if (db[i] >= 1) {
        ai = Norm(1, t2[i]);

	ci = 1 / ai + 1;
	di = t2[i] / ai + 1;
      } 
      else {
        bi = x[i] - u[i];
        ai = Norm(bi, f[i]);
	ai = PetscMax(TAO_EPSILON, ai);

	ci = bi / ai + 1;
	di = f[i] / ai + 1;
      }

      if (da[i] >= 1) {
        bi = ci + di*t2[i];
        ai = Norm(1, bi);
	
	bi = bi / ai - 1;
	ai = 1 / ai - 1;
      } 
      else {
	ei = Fischer(u[i] - x[i], -f[i]);
	ai = Norm(x[i] - l[i], ei);
	ai = PetscMax(TAO_EPSILON, ai);

	bi = ei / ai - 1;
	ai = (x[i] - l[i]) / ai - 1;
      }

      da[i] = ai + bi*ci;
      db[i] = bi*di;
    }
  }

  info = VecRestoreArray(DA, &da); CHKERRQ(info);
  info = VecRestoreArray(DB, &db); CHKERRQ(info);

  info = VecRestoreArray(X, &x); CHKERRQ(info);
  info = VecRestoreArray(F, &f); CHKERRQ(info);
  info = VecRestoreArray(L, &l); CHKERRQ(info);
  info = VecRestoreArray(U, &u); CHKERRQ(info);
  info = VecRestoreArray(T2, &t2); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SFischer"
inline static PetscScalar SFischer(PetscScalar a, PetscScalar b, PetscScalar c)
{

#ifdef TODD

  if (PetscAbsScalar(a) > PetscAbsScalar(b)) {
    return sqrt(a*a + b*b + 2.0*c*c) - a - b;
  }
  return sqrt(a*a + b*b + 2.0*c*c) - b - a;

#else

   // Method suggested by Bob Vanderbei

   if (a + b <= 0) {
     return sqrt(a*a + b*b + 2.0*c*c) - (a + b);
   }
   return 2.0*(c*c - a*b) / (sqrt(a*a + b*b + 2.0*c*c) + (a + b));

#endif

}

#undef __FUNCT__
#define __FUNCT__ "SNorm"
inline static PetscScalar SNorm(PetscScalar a, PetscScalar b, PetscScalar c)
{
  return sqrt(a*a + b*b + 2.0*c*c);
}

#undef __FUNCT__
#define __FUNCT__ "MatD_SFischer"
int MatD_SFischer(Mat m, Vec X, Vec F, Vec L, Vec U, double mu,
                  Vec T1, Vec T2, Vec DA, Vec DB, Vec DM)
{
  PetscScalar *x, *f, *l, *u, *da, *db, *dm;
  PetscScalar ai, bi, ci, di, ei, fi;
  int info, i;
  PetscInt n, low[7], high[7];

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X, VEC_COOKIE,2);
  PetscValidHeaderSpecific(F, VEC_COOKIE,3);
  PetscValidHeaderSpecific(L, VEC_COOKIE,4);
  PetscValidHeaderSpecific(U, VEC_COOKIE,5);
  PetscValidHeaderSpecific(DA, VEC_COOKIE,9);
  PetscValidHeaderSpecific(DB, VEC_COOKIE,10);
  PetscValidHeaderSpecific(DM, VEC_COOKIE,11);

  info = VecGetOwnershipRange(X, low, high); CHKERRQ(info);
  info = VecGetOwnershipRange(F, low + 1, high + 1); CHKERRQ(info);
  info = VecGetOwnershipRange(L, low + 2, high + 2); CHKERRQ(info);
  info = VecGetOwnershipRange(U, low + 3, high + 3); CHKERRQ(info);
  info = VecGetOwnershipRange(DA, low + 4, high + 4); CHKERRQ(info);
  info = VecGetOwnershipRange(DB, low + 5, high + 5); CHKERRQ(info);
  info = VecGetOwnershipRange(DM, low + 6, high + 6); CHKERRQ(info);

  for (i = 1; i < 7; i++) {
    if (low[0] != low[i] || high[0] != high[i])
      SETERRQ(1,"Vectors must be identically loaded over processors");
  }

  info = VecGetArray(X, &x); CHKERRQ(info);
  info = VecGetArray(F, &f); CHKERRQ(info);
  info = VecGetArray(L, &l); CHKERRQ(info);
  info = VecGetArray(U, &u); CHKERRQ(info);
  info = VecGetArray(DA, &da); CHKERRQ(info);
  info = VecGetArray(DB, &db); CHKERRQ(info);
  info = VecGetArray(DM, &dm); CHKERRQ(info);

  info = VecGetLocalSize(X, &n); CHKERRQ(info);

  for (i = 0; i < n; i++) {
    if ((l[i] <= -TAO_INFINITY) && (u[i] >= TAO_INFINITY)) {
      da[i] = -mu;
      db[i] = -1;
      dm[i] = -x[i];
    } 
    else if (l[i] <= -TAO_INFINITY) {
      bi = u[i] - x[i];
      ai = SNorm(bi, f[i], mu);
      ai = PetscMax(TAO_EPSILON, ai);

      da[i] = bi / ai - 1;
      db[i] = -f[i] / ai - 1;
      dm[i] = 2.0 * mu / ai;
    } 
    else if (u[i] >=  TAO_INFINITY) {
      bi = x[i] - l[i];
      ai = SNorm(bi, f[i], mu);
      ai = PetscMax(TAO_EPSILON, ai);

      da[i] = bi / ai - 1;
      db[i] = f[i] / ai - 1;
      dm[i] = 2.0 * mu / ai;
    } 
    else if (l[i] == u[i]) {
      da[i] = -1;
      db[i] = 0;
      dm[i] = 0;
    } 
    else {
      bi = x[i] - u[i];
      ai = SNorm(bi, f[i], mu);
      ai = PetscMax(TAO_EPSILON, ai);

      ci = bi / ai + 1;
      di = f[i] / ai + 1;
      fi = 2.0 * mu / ai;

      ei = SFischer(u[i] - x[i], -f[i], mu);
      ai = SNorm(x[i] - l[i], ei, mu);
      ai = PetscMax(TAO_EPSILON, ai);

      bi = ei / ai - 1;
      ai = (x[i] - l[i]) / ai - 1;

      da[i] = ai + bi*ci;
      db[i] = bi*di;
      dm[i] = ei + bi*fi;
    }
  }

  info = VecRestoreArray(X, &x); CHKERRQ(info);
  info = VecRestoreArray(F, &f); CHKERRQ(info);
  info = VecRestoreArray(L, &l); CHKERRQ(info);
  info = VecRestoreArray(U, &u); CHKERRQ(info);
  info = VecRestoreArray(DA, &da); CHKERRQ(info);
  info = VecRestoreArray(DB, &db); CHKERRQ(info);
  info = VecRestoreArray(DM, &dm); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#endif
