#ifndef PETSCDAAPPHEADER_H
#define PETSCDAAPPHEADER_H


/* 
   This file declares the interface to optimization application defined
   over a set of PETSc Distributed Array (DA) objects
*/
#include "tao_general.h"
#include "taoapp.h"
#include "petscksp.h"
#include "petscda.h"

typedef struct _p_TAOAPPLICATION* TAO_APPLICATION;

extern int TaoAppSetDAApp(TAO_APPLICATION, DA*, PetscInt);
extern int DAApplicationCreate(MPI_Comm, DA*, PetscInt, TAO_APPLICATION*);
extern int TaoAppDestroyDAApp(TAO_APPLICATION);
extern int DAAppSetHessianRoutine(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,Vec,Mat,void*),void *);
extern int DAAppSetObjectiveAndGradientRoutine(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,Vec,double*,Vec, void*),void *);
extern int DAAppSetGradientRoutine(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,Vec,Vec, void*),void *);
extern int DAAppSetObjectiveRoutine(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,Vec,double*, void*),void *);
extern int DAAppSetVariableBoundsRoutine(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,Vec,Vec, void*),void *);
extern int DAAppSetConstraintRoutine(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,Vec,Vec, void*),void *);
extern int DAAppSetJacobianRoutine(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,Vec,Mat,void*),void *);

extern int DAAppSetBeforeMonitor(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,PetscInt,void*),void *);
extern int DAAppSetAfterMonitor(TAO_APPLICATION, int (*)(TAO_APPLICATION,DA,PetscInt,void*),void *);
extern int DAAppSetInitialSolution(TAO_APPLICATION,Vec);
extern int DAAppSetHessianMat(TAO_APPLICATION);
extern int DAAppSetMatType(TAO_APPLICATION, const MatType);
extern int DAAppSetQuadraticObjective(TAO_APPLICATION, PetscTruth);
extern int DAAppSetMatStructure(TAO_APPLICATION, MatStructure);
extern int DAAppGetSolution(TAO_APPLICATION, PetscInt, Vec *);
extern int DAAppGetHessianMat(TAO_APPLICATION, PetscInt, Mat *);
extern int DAAppGetVariableBounds(TAO_APPLICATION, PetscInt, Vec *, Vec*);
extern int DAAppGetInterpolationMatrix(TAO_APPLICATION, PetscInt, Mat *, Vec *);
extern int DAAppGetDA(TAO_APPLICATION, PetscInt, DA *);
extern int DAAppGetCurrentLevel(TAO_APPLICATION, PetscInt*);
extern int DAAppGetNumberOfDAGrids(TAO_APPLICATION, PetscInt*);

extern int DAAppSetOptions(TAO_APPLICATION);

/* For developers to access the internal data structure */
typedef struct _p_DA_APPLICATION* DA_APPLICATION;
extern int TaoAppDAApp(TAO_APPLICATION, DA_APPLICATION *);

/* Define objective on a finite element -- Additional support for DAApp, but separate from previous routines */
extern int DAAppSetElementObjectiveAndGradientRoutine(TAO_APPLICATION, int (*)(PetscInt[2],double[4],double*,double[4],void*), PetscInt, void *);
extern int DAAppSetElementHessianRoutine(TAO_APPLICATION, int (*)(PetscInt[2],double[4],double[4][4],void*), PetscInt, void *);

/* Use multigrid preconditioning -- Additional support, but separate from previous routines */
extern int DAAppUseMultigrid(TAO_APPLICATION,PetscInt);
extern int DAAppSetupMultigrid(TAO_APPLICATION, PetscInt);
extern int DAAppSetMultiGridOptions(TAO_APPLICATION);

extern int DAAppPrintStageTimes(TAO_APPLICATION);
extern int DAAppTimeMonitorAfter(TAO_APPLICATION,DA,PetscInt,void*);

extern int DAAppInterpolationMonitor(TAO_APPLICATION,DA,PetscInt,void*);
extern int DAAppPrintInterpolationError(TAO_APPLICATION);

typedef struct _p_TAO_SOLVER* TAO_SOLVER;
extern int TaoDAAppSolve(TAO_APPLICATION, TAO_SOLVER);

/* Use ADIC for objetive function -- Additional but separate  */
#include "src/petsctao/gridapplication/daadic/taodaadic.h"

#endif
