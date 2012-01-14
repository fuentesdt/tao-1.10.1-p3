#include "tao_general.h"
#include "taodaapplication.h"     /* taodaapplication.h */
#include "src/petsctao/include/taopetsc.h"

#include "../../application/taoapp/taoapp_petsc.h"
#include "../interface/daapp_impl.h"
#include "multigridmat.h"

#include "petscmg.h"

class TaoMultiGridApplication: public TaoPetscApplication {

 protected:

 public:

  TaoMultiGridApplication(MPI_Comm);
  TaoMultiGridApplication();

  PetscInt coarselevel;
  int EvaluateHessian(TaoVec *xx, TaoMat *HH);

};


#undef __FUNCT__  
#define __FUNCT__ "TaoMultiGridApplication::TaoMultiGridApplication"
TaoMultiGridApplication :: TaoMultiGridApplication(MPI_Comm mpicomm) : TaoPetscApplication(mpicomm) {
  Mat M=0;
  this->coarselevel=0;
  this->taoh=new TaoMultiGridMatPetsc(M);
  this->taoj=new TaoMultiGridMatPetsc(M);
  return;
}

#undef __FUNCT__  
#define __FUNCT__ "TaoMultiGridApplication::TaoMultiGridApplication"
TaoMultiGridApplication::TaoMultiGridApplication() : TaoPetscApplication(PETSC_COMM_WORLD){
  
  Mat M=0;
  this->coarselevel=0;
  this->taoh=new TaoMultiGridMatPetsc(M);
  this->taoj=new TaoMultiGridMatPetsc(M);
  return;
}


#undef __FUNCT__  
#define __FUNCT__ "TaoMultiGridApplication::EvaluateHessian"
int TaoMultiGridApplication::EvaluateHessian(TaoVec *xx, TaoMat *HH)
{
  TaoVecPetsc *px = dynamic_cast <TaoVecPetsc *> (xx);
  TaoMatPetsc *pH = dynamic_cast <TaoMatPetsc *> (HH);

  int     info;
  PetscInt i,clevel,cclevel;
  Vec X = px->GetVec();
  Mat H, HPre;
  MatStructure pflag;
  DA_APPLICATION daapp;

  PetscValidHeaderSpecific(papp,TAO_APP_COOKIE,0);
  info=TaoAppDAApp(this->papp,&daapp); CHKERRQ(info);

  PetscStackPush("TAO User DA Hessian");
  info=pH->GetMatrix(&H,&HPre,&pflag);CHKERRQ(info);
  info=DAAppGetCurrentLevel(this->papp,&clevel); CHKERRQ(info);
  cclevel=PetscMin(clevel,this->coarselevel);
  cclevel=PetscMax(cclevel,0);

  for (i=clevel; i>=cclevel; i--){
    daapp->currentlevel=i;
    info=TaoAppSetColoring(this->papp,daapp->grid[i].coloring); CHKERRQ(info);
    info = PetscInfo1(daapp,"TAO hessian evaluation at DA_APPLICATION object, level %d.\n",i); CHKERRQ(info);
    info = TaoAppComputeHessian(this->papp, X, &H, &HPre, &pflag); CHKERRQ(info);
    if (i==clevel){  
      info = pH->SetMatrix(H,HPre,pflag);CHKERRQ(info);
    }

    if (i>cclevel){
      info=MatMultTranspose(daapp->grid[i].Interpolate,X,daapp->grid[i-1].R);CHKERRQ(info);
      info=VecPointwiseMult(daapp->grid[i-1].R,daapp->grid[i].CScale,daapp->grid[i-1].R);CHKERRQ(info);
      X=daapp->grid[i-1].R;
      H=daapp->grid[i-1].H;
      HPre=H;
    }
  }
  daapp->currentlevel=clevel;
  info=TaoAppSetColoring(this->papp,daapp->grid[clevel].coloring); CHKERRQ(info);
  PetscStackPop;

  PetscFunctionReturn(0);
}

int DAAppSetMultigridKSP(GridCtx *, PetscInt, KSP);
int DAAppSetupMultigridMonitor(TAO_APPLICATION, DA, PetscInt, void *);
int TaoAppGetMultiGridApplication(TAO_APPLICATION, TaoMultiGridApplication *);
static char TaoPetscPCMGDAAppXXX[] = "TaoPetscPCMGDAApp";
static char TaoPetscAppXXX[] = "TaoPetscApp";


#undef __FUNCT__  
#define __FUNCT__ "DAAppDestroyTaoAppXXX"
static int DAAppDestroyTaoAppXXX(void*ctx){
  TaoMultiGridApplication *mctx=(TaoMultiGridApplication*)ctx;
  PetscFunctionBegin;
  if (mctx){ delete mctx;}
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetMultiGridApplication"
int TaoAppGetMultiGridApplication(TAO_APPLICATION daapplication, TaoMultiGridApplication **mgdaapp){
  int info;
  PetscInt clevel,ii;
  Mat H,HH;
  MatStructure         pflag=SAME_NONZERO_PATTERN;
  DA_APPLICATION       daapp;
  TaoMultiGridApplication*   daapp2;
  TaoPetscApplication  *pppappp;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(daapplication,TAO_APP_COOKIE,1);
  info=TaoAppQueryForObject(daapplication,TaoPetscPCMGDAAppXXX,(void**)&daapp2); CHKERRQ(info);
  if (daapp2==0){
    daapp2=new TaoMultiGridApplication(PETSC_COMM_WORLD);
    daapp2->papp=daapplication;
    info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);

    info=TaoAppQueryForObject(daapplication,TaoPetscAppXXX,(void**)&pppappp); CHKERRQ(info);
    if (pppappp && pppappp->taoh && 0==1){
      info=pppappp->taoh->GetMatrix(&H,&HH,&pflag); CHKERRQ(info);
      info=daapp2->taoh->SetMatrix(H,HH,pflag);CHKERRQ(info);
    } else {
      info=DAAppGetCurrentLevel(daapplication,&clevel); CHKERRQ(info);
      info=daapp2->taoh->SetMatrix(daapp->grid[clevel].H,daapp->grid[clevel].H,pflag);CHKERRQ(info);
    }
    if (pppappp){
      daapp2->tao=pppappp->tao;
    }
    info=TaoAppAddObject(daapplication,TaoPetscPCMGDAAppXXX,(void*)daapp2,&ii); CHKERRQ(info);
    info=TaoAppQueryRemoveObject(daapplication,TaoPetscAppXXX); CHKERRQ(info);
    info=TaoAppAddObject(daapplication,TaoPetscAppXXX,(void*)daapp2,&ii); CHKERRQ(info);
    info=TaoAppSetDestroyRoutine(daapplication,DAAppDestroyTaoAppXXX, (void*)daapp2); CHKERRQ(info);
  }
  *mgdaapp=daapp2;
  PetscFunctionReturn(0);
}

static int UsePCMGGG=0;
#undef __FUNCT__
#define __FUNCT__ "DAAppUseMultigrid"
/*@
  DAAppUseMultigrid - Set the preconditioner for the linear solver to be an algebraic multigrid.

  Collective on TAO_APPLICATION
  
  Input Parameters:
+  daapplication - the DA Application object
-  coarselevel - the coarsest grid to be used in the multigrid preconditioner. (Grid 0 is the coarsest grid.

   Level: intermediate

   Options Database Key:
+  -tao_da_multigrid - use multigrid linear solver
-  -ksp_view - view the linear solver

Note:
  This function should be called after DAAppSetHessianRoutine();

Note:
  This function should be called before each optimization solver as part of the DAAppMonitor

Note:
  Multigrid functionality is still under developement for good performance.

.seealso: TaoAppGetKSP(), DAAppSetupMultigrid()

   Options Database Key:
.  -tao_da_multigrid

.keywords: Linear Solver, Multigrid, DA, KSP

@*/
int DAAppUseMultigrid(TAO_APPLICATION daapplication, PetscInt coarselevel){
  int info;
  PetscFunctionBegin;
  UsePCMGGG=coarselevel;
  info=DAAppSetBeforeMonitor(daapplication,DAAppSetupMultigridMonitor,0); 
  CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DAAppSetupMultigridMonitor"
int DAAppSetupMultigridMonitor(TAO_APPLICATION daapplication, DA da, PetscInt level, void *ctx){
  int info;
  PetscFunctionBegin;
  info = DAAppSetupMultigrid(daapplication,UsePCMGGG); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DAAppSetupMultigrid"
/*@
   DAAppSetupMultigrid - Sets up the multigrid preconditioner.

   Collective on TAO_APPLICATION

   Input Parameters:
+  daapplication - the TAO_APPLICATION context
-  coarselevel - the coarsest grid that should be used in the multigrid (>=0)

   Note:
   Usually the coarselevel is set to 0;

   Level: intermediate

.seealso:  DAAppUseMultigrid(), TaoAppSetMonitor()
@*/
int DAAppSetupMultigrid(TAO_APPLICATION daapplication, PetscInt coarselevel){
  int info;
  PetscInt nn,level;
  TaoMultiGridMatPetsc *MMM;
  DA_APPLICATION daapp;
  TaoMultiGridApplication *mgdaapp;
  TAO_SOLVER tao;
  KSP  ksp;
 
  PetscFunctionBegin;
  info=DAAppGetCurrentLevel(daapplication,&level); CHKERRQ(info);
  info=TaoAppDAApp(daapplication,&daapp); CHKERRQ(info);
  if (daapp->grid[level].mgrid==0 &&  coarselevel<level){
    info=TaoAppGetMultiGridApplication(daapplication,&mgdaapp); CHKERRQ(info);
    mgdaapp->coarselevel=PetscMax(0,coarselevel);
    nn=level-coarselevel+1;
    if ( mgdaapp->taoh ){
      MMM=(TaoMultiGridMatPetsc*)(mgdaapp->taoh);
      info=MMM->SetUpMultiGrid(daapp->grid+coarselevel,nn); CHKERRQ(info);
    } else if ( mgdaapp->taoj ){
      MMM=(TaoMultiGridMatPetsc*)(mgdaapp->taoj);
      info=MMM->SetUpMultiGrid(daapp->grid+coarselevel,nn); CHKERRQ(info);
      mgdaapp->taoj=MMM;
    }

    info=TaoAppGetTaoSolver(daapplication,&tao); CHKERRQ(info);
    info=TaoSetupApplicationSolver(daapplication,tao); CHKERRQ(info);
    info=TaoAppGetKSP(daapplication,&ksp); CHKERRQ(info);
    if (ksp){
      PetscValidHeaderSpecific(ksp,KSP_COOKIE,1);
      info = DAAppSetMultigridKSP(daapp->grid+coarselevel,nn,ksp);CHKERRQ(info);
    }
    daapp->grid[level].mgrid=1;
  }
  PetscFunctionReturn(0);
}


  
#undef __FUNCT__
#define __FUNCT__ "DAAppSetMultigridKSP"
int DAAppSetMultigridKSP(GridCtx* dagrid, PetscInt nlevels, KSP ksp){
  int info;
  PetscInt i;
  PC pc,pc2;
  KSP sksp;
  const KSPType ksptype;
  const PCType pctype;
  PCMGType form;

  PetscFunctionBegin;
  if (ksp==0 || nlevels<=1){
    PetscFunctionReturn(0);
  }

  PetscValidHeaderSpecific(ksp,KSP_COOKIE,3);
  info = PetscInfo(ksp,"Set multigrid precondition in optimization solver.\n"); CHKERRQ(info);
  
  info = KSPGetType(ksp,&ksptype); CHKERRQ(info);
  
  info = KSPGetPC(ksp,&pc); CHKERRQ(info);
  info = PCGetType(pc,&pctype); CHKERRQ(info);
  info = PCSetType(pc,PCMG); CHKERRQ(info);
  
  info = PCMGSetLevels(pc,nlevels,PETSC_NULL); CHKERRQ(info);

  for (i=0; i<nlevels; i++) {
    
    info = PCMGGetSmoother(pc,i,&sksp); CHKERRQ(info);
    info = KSPSetType(sksp,ksptype); CHKERRQ(info);

    info = KSPGetPC(sksp,&pc2); CHKERRQ(info);
    info = PCSetType(pc2,PCJACOBI); CHKERRQ(info);

    form=PC_MG_MULTIPLICATIVE;
    info=PCMGSetType(pc,form);CHKERRQ(info);

    info = KSPSetOperators(sksp,dagrid[i].H,dagrid[i].H,SAME_NONZERO_PATTERN); CHKERRQ(info);
    
    if (i<nlevels-1){
      info = PCMGSetX(pc,i,dagrid[i].X); CHKERRQ(info);
      info = PCMGSetRhs(pc,i,dagrid[i].RHS); CHKERRQ(info);      
    }
    if (i>0){
      info = PCMGSetR(pc,i,dagrid[i].R); CHKERRQ(info);
    }
    info = PCMGSetResidual(pc,i,PCMGDefaultResidual,dagrid[i].H); CHKERRQ(info);
  }
  for (i=1; i<nlevels; i++) {
    info = PCMGSetInterpolation(pc,i,dagrid[i].Interpolate); CHKERRQ(info);
    info = PCMGSetRestriction(pc,i,dagrid[i].Interpolate); CHKERRQ(info);
  }
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAAppSetQuadraticObjective"
/* @
  DAAppSetQuadraticObjective - identify the objective function as quadratic or not.

  Collective on TAO_APPLICATION

  Input Parameters:
+  daapplication - the DA Application object
-  flag - indicates whether the objective is quadratic or not

   Level: intermediate

   Options Database Key:
.  -tao_da_quadratic - use multigrid linear solver

   Note:
   The default value is PETSC_FALSE.

   Note:
   If quadratic, consider setting the flag used in KSP to SAME_PRECONDITIONER

   Note:
   If quadratic, this routine reduces the number of Hessian evaluations done when using
   the multigrid preconditioner.

.seealso: DAAppSetMatStructure(),  DAAppUseMultigrid()


.keywords: DA, Objective Function

@ */
int DAAppSetQuadraticObjective(TAO_APPLICATION daapplication, PetscTruth pflag){
  int info;
  TaoMultiGridApplication *mgdaapp;

  PetscFunctionBegin;
  info=TaoAppGetMultiGridApplication(daapplication,&mgdaapp); CHKERRQ(info);
  //  mgdaapp->IsQuadratic=pflag;
  if (pflag==PETSC_TRUE){
    info = PetscInfo(daapplication,"User labels the objective function as quadratic.\n"); CHKERRQ(info);
    info = PetscInfo(daapplication,"Set KSP MatStructure to SAME_PRECONDITIONER.\n"); CHKERRQ(info);
    info=DAAppSetMatStructure(daapplication, SAME_PRECONDITIONER);CHKERRQ(info);
  } else {
    info = PetscInfo(daapplication,"User labels the objective function as NOT quadratic.\n"); CHKERRQ(info);
    info = PetscInfo(daapplication,"Set KSP MatStructure to SAME_NONZERO_PATTERN.\n"); CHKERRQ(info);
    info=DAAppSetMatStructure(daapplication, SAME_NONZERO_PATTERN);CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "DAAppSetMultiGridOptions"
/* @
  DAAppSetMultiGridOptions - Sets various multigrid options to be used in this application
  and the TAO solver.

   Collective on TAO_APPLICATION

   Input Parameters:
.  daapplication - the DA Application object

   Level: beginner

.keywords:  options

.seealso: TaoDAAppSetOptions();

@ */
int DAAppSetMultiGridOptions(TAO_APPLICATION daapplication){
  int info;
  PetscInt nlevel;

  PetscTruth flg1=PETSC_FALSE;

  PetscFunctionBegin;
  info = PetscInfo(daapplication,"TaoDAAppSetMultiGridOptions(): Reading command line for options\n"); CHKERRQ(info);

  flg1=PETSC_FALSE;
  info = PetscOptionsInt("-tao_da_multigrid","use multigrid linear solver","DAAppUseMultigrid",
			 PETSC_FALSE,&nlevel,&flg1);CHKERRQ(info);
  if (flg1) {
    info=DAAppUseMultigrid(daapplication,nlevel);CHKERRQ(info);
    info = PetscInfo(daapplication,"TaoDAAppSetOptions: Use Multigrid linear solver \n"); CHKERRQ(info);
  }
  /*
  flg1=PETSC_FALSE,flg2=PETSC_FALSE;
  info = PetscOptionsTruth("-tao_da_quadratic","Identify the objective function as quadratic",
			     "DAAppSetQuadraticObjective",PETSC_FALSE,&flg2,&flg1);CHKERRQ(info);
  if (flg1) {
    info = DAAppSetQuadraticObjective(daapplication, flg2);CHKERRQ(info); 
  }
  */

  PetscFunctionReturn(0);
}

