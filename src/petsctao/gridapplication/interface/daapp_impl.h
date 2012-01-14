#ifndef PETSCDAAPPIMPL_H
#define PETSCDAAPPIMPL_H

#include "petscmat.h"
#include "petscda.h"
#include "dagridctx.h"

typedef struct _p_DA_APPLICATION* DA_APPLICATION;
typedef struct _p_TAOAPPLICATION* TAO_APPLICATION;


#define PETSCDAAPPMAXGRIDS 20
#define MAX_DAAP_MONITORS 10

struct _p_DA_APPLICATION {

  PETSCHEADER(int);

  GridCtx grid[PETSCDAAPPMAXGRIDS];

  PetscInt ndamax;        /* Max number of levels */
  PetscInt nda;           /* Number of levels in current application */
  PetscInt currentlevel;  /* The current level being solved */

  PetscTruth IsComplementarity;
  MatStructure kspflag;
  char HessianMatrixType[20];
  
  /* Function Gradient Evaluation over entire DA */
  int  (*computedafunction)(TAO_APPLICATION,DA,Vec,double*,void*); 
  int  (*computedagradient)(TAO_APPLICATION,DA,Vec,Vec,void*); 
  int  (*computedafunctiongradient)(TAO_APPLICATION,DA,Vec,double*,Vec,void*); 
  void *usrdafctx;
  void *usrdagctx;
  void *usrdafgctx;

  /* Hessian Evaluation over entire DA */
  int  (*computedahessian)(TAO_APPLICATION,DA,Vec,Mat,void*);
  void *usrdahctx;

  /* Evaluate bounds over entire DA */
  void* bounddactx;
  int  (*computedabounds)(TAO_APPLICATION,DA,Vec,Vec,void*);

  /* User monitors before and after each Optimization Solve */
  PetscInt   nbeforemonitors;
  void *beforemonitorctx[MAX_DAAP_MONITORS];
  int  (*beforemonitor[MAX_DAAP_MONITORS])(TAO_APPLICATION,DA,PetscInt,void*); 

  PetscInt   naftermonitors;
  void *aftermonitorctx[MAX_DAAP_MONITORS];
  int  (*aftermonitor[MAX_DAAP_MONITORS])(TAO_APPLICATION,DA,PetscInt,void*); 

};

//#include "taodaapplication.h"

#endif


