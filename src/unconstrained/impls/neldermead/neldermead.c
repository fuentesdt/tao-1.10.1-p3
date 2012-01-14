#include "neldermead.h"
#ifdef TAO_USE_PETSC
#include "src/petsctao/vector/taovec_petsc.h"
#endif

int NelderMeadSort(TAO_NelderMead *nm);
int NelderMeadReplace(TAO_NelderMead *nm, int index, TaoVec *Xmu, double f);
/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetUp_NelderMead"
int TaoSetUp_NelderMead(TAO_SOLVER tao, void *solver)
{
  TAO_NelderMead *nm = (TAO_NelderMead *)solver;
  int info;
  TaoInt size;
  TaoVec *X;

  TaoFunctionBegin;

  
  info = TaoGetSolution(tao,&X); CHKERRQ(info);
  info = X->GetDimension(&size); CHKERRQ(info);
  nm->N = size;
  nm->oneOverN = 1.0/size;
  info = X->CloneVecs(nm->N+1,&nm->simplex);
  nm->f_values = new double[nm->N+1];
  nm->indices =new int[nm->N+1];
  info = X->Clone(&nm->Xbar); CHKERRQ(info);
  info = X->Clone(&nm->Xmur); CHKERRQ(info);
  info = X->Clone(&nm->Xmue); CHKERRQ(info);
  info = X->Clone(&nm->Xmuc); CHKERRQ(info);
  info = TaoSetLagrangianGradientVector(tao,0);CHKERRQ(info);
  info = TaoSetStepDirectionVector(tao,0);CHKERRQ(info);

  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoSetDown_NelderMead"
int TaoSetDown_NelderMead(TAO_SOLVER tao, void *solver)
{
  TAO_NelderMead *nm = (TAO_NelderMead*)solver;
  int info;
  int i;
  TaoFunctionBegin;
  for (i=0;i<nm->N+1;i++) {
    TaoVecDestroy( nm->simplex[i]);
  }
  delete [] nm->simplex;

  info = TaoVecDestroy(nm->Xmuc); CHKERRQ(info);
  info = TaoVecDestroy(nm->Xmue); CHKERRQ(info);
  info = TaoVecDestroy(nm->Xmur); CHKERRQ(info);
  info = TaoVecDestroy(nm->Xbar); CHKERRQ(info);
  
  delete [] nm->indices;
  delete [] nm->f_values;
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSetOptions_NelderMead"
int TaoSetOptions_NelderMead(TAO_SOLVER tao, void *solver)
{
  
  TAO_NelderMead *nm = (TAO_NelderMead*)solver;
  TaoTruth flg;
  int info;
  
  TaoFunctionBegin;
  info = TaoOptionsHead("Nelder-Mead options"); CHKERRQ(info);

  info = TaoOptionDouble("-tao_nm_lamda","initial step length","",nm->lamda,&nm->lamda,&flg);

  CHKERRQ(info);

  info = TaoOptionDouble("-tao_nm_mu","mu","",nm->mu_oc,&nm->mu_oc,&flg); CHKERRQ(info);
  nm->mu_ic = -nm->mu_oc;
  nm->mu_r = nm->mu_oc*2.0;
  nm->mu_e = nm->mu_oc*4.0;

  info = TaoOptionsTail(); CHKERRQ(info);
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoView_NelderMead"
int TaoView_NelderMead(TAO_SOLVER tao, void *solver)
{
  //int info;

  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "TaoSolve_NelderMead"
int TaoSolve_NelderMead(TAO_SOLVER tao, void *solver)
{
  int info;
  TAO_NelderMead *nm = (TAO_NelderMead*)solver;
  TaoTerminateReason reason;
  TaoVec *xx;
  double *x;
  double step=0.0;
  int iter=0,i;
  TaoVec *Xmur=nm->Xmur, *Xmue=nm->Xmue, *Xmuc=nm->Xmuc, *Xbar=nm->Xbar;
  double fr,fe,fc;
  int shrink;
  
  
  TaoFunctionBegin;
  info = TaoGetSolution(tao,&xx); CHKERRQ(info);
  info = nm->simplex[0]->CopyFrom(xx); CHKERRQ(info);
  info = TaoComputeMeritFunction(tao,nm->simplex[0],&nm->f_values[0]); CHKERRQ(info);
  nm->indices[0] = 0;
  for (i=1;i<nm->N+1;i++){
    info = nm->simplex[i]->CopyFrom(xx); CHKERRQ(info);
#ifdef TAO_USE_PETSC
    PetscInt low,high;
    TaoVecPetsc *tvxx  = dynamic_cast<TaoVecPetsc*>(nm->simplex[i]);
    if (!tvxx) {SETERRQ(1,"Could not cast TaoVec to TaoVecPetsc");}
    Vec px = tvxx->GetVec();
    info = VecGetOwnershipRange(px,&low,&high); CHKERRQ(info);
    if (i-1 >= low && i-1 < high) {
      info = VecGetArray(px,&x); CHKERRQ(info);
      x[i-1-low] += nm->lamda;
      info = VecRestoreArray(px,&x); CHKERRQ(info);
    }


#else
    TaoInt dim;
    if (i>0) {
      
      info = nm->simplex[i]->GetArray(&x,&dim); CHKERRQ(info);
      x[i-1] += nm->lamda;
      info = nm->simplex[i]->RestoreArray(&x,&dim); CHKERRQ(info);
    }
#endif
    info = TaoComputeMeritFunction(tao,nm->simplex[i],&nm->f_values[i]); CHKERRQ(info);
    nm->indices[i] = i;
  }

  // Xbar  = (Sum of all simplex vectors - worst vector)/N
  info = NelderMeadSort(nm); CHKERRQ(info);
  info = Xbar->SetToZero(); CHKERRQ(info);
  for (i=0;i<nm->N;i++) {
    info = Xbar->Axpy(1.0,nm->simplex[nm->indices[i]]);
  }
  info = Xbar->Scale(nm->oneOverN);

  reason = TAO_CONTINUE_ITERATING;
  while (1) {
    shrink = 0;

    info = xx->CopyFrom(nm->simplex[nm->indices[0]]); CHKERRQ(info);
    info = TaoMonitor(tao,iter++,nm->f_values[nm->indices[0]],nm->f_values[nm->indices[nm->N]]-nm->f_values[nm->indices[0]],0,step,&reason); CHKERRQ(info);
    if (reason != TAO_CONTINUE_ITERATING) break;

    
    
    //x(mu) = (1 + mu)Xbar - mu*X_N+1
    info = Xmur->Waxpby(1+nm->mu_r,Xbar,-nm->mu_r,
		       nm->simplex[nm->indices[nm->N]]); CHKERRQ(info);
    info = TaoComputeMeritFunction(tao,Xmur,&fr); CHKERRQ(info);


    if (nm->f_values[nm->indices[0]] <= fr && fr < nm->f_values[nm->indices[nm->N-1]]) {
      // reflect
      info = PetscInfo(0,"Reflect\n"); CHKERRQ(info);
      info = NelderMeadReplace(nm,nm->indices[nm->N],Xmur,fr); CHKERRQ(info);
    }

    else if (fr < nm->f_values[nm->indices[0]]) {
      // expand
      info = PetscInfo(0,"Expand\n"); CHKERRQ(info);
      info = Xmue->Waxpby(1+nm->mu_e,Xbar,-nm->mu_e,
			  nm->simplex[nm->indices[nm->N]]); CHKERRQ(info);
      info = TaoComputeMeritFunction(tao,Xmue,&fe); CHKERRQ(info);
      if (fe < fr) {
	info = NelderMeadReplace(nm,nm->indices[nm->N],Xmue,fe); CHKERRQ(info);
      } else {
	info = NelderMeadReplace(nm,nm->indices[nm->N],Xmur,fr); CHKERRQ(info);
      }

    } else if (nm->f_values[nm->indices[nm->N-1]] <= fr && fr < nm->f_values[nm->indices[nm->N]]) {
      //outside contraction
      info = PetscInfo(0,"Outside Contraction\n"); CHKERRQ(info);
      info = Xmuc->Waxpby(1+nm->mu_oc,Xbar,-nm->mu_oc,
			  nm->simplex[nm->indices[nm->N]]); CHKERRQ(info);
      info = TaoComputeMeritFunction(tao,Xmuc,&fc); CHKERRQ(info);
      if (fc <= fr) {
	info = NelderMeadReplace(nm,nm->indices[nm->N],Xmuc,fc); CHKERRQ(info);
      }	else 
	shrink=1;
    } else {
      //inside contraction
      info = PetscInfo(0,"Inside Contraction\n"); CHKERRQ(info);
      info = Xmuc->Waxpby(1+nm->mu_ic,Xbar,-nm->mu_ic,
			  nm->simplex[nm->indices[nm->N]]); CHKERRQ(info);
      info = TaoComputeMeritFunction(tao,Xmuc,&fc); CHKERRQ(info);
      if (fc < nm->f_values[nm->indices[nm->N]]) {
	info = NelderMeadReplace(nm,nm->indices[nm->N],Xmuc,fc); CHKERRQ(info);
      } else
	shrink = 1;
    }

    if (shrink) {
      info = PetscInfo(0,"Shrink\n"); CHKERRQ(info);
      printf("shrink!\n");
      
      for (i=1;i<nm->N+1;i++) {
	info = nm->simplex[nm->indices[i]]->Axpby(1.5,nm->simplex[nm->indices[0]],-.5); CHKERRQ(info);
	info = TaoComputeMeritFunction(tao,nm->simplex[nm->indices[i]],
				 &nm->f_values[nm->indices[i]]); CHKERRQ(info);
      }

      info = Xbar->Axpby(1.5*nm->oneOverN,nm->simplex[nm->indices[0]],-0.5); CHKERRQ(info);
      // Add last vector's fraction of average
      info = nm->Xbar->Axpy(nm->oneOverN,
			    nm->simplex[nm->indices[nm->N]]); CHKERRQ(info);
      info = NelderMeadSort(nm);
      // Subtract new last vector from average
      info = nm->Xbar->Axpy(-nm->oneOverN,
			    nm->simplex[nm->indices[nm->N]]); CHKERRQ(info);
    }
    
    
  }
  TaoFunctionReturn(0);
}

/* ---------------------------------------------------------- */
EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "TaoCreate_NelderMead"
int TaoCreate_NelderMead(TAO_SOLVER tao)
{
  TAO_NelderMead *nm;
  int info;

  TaoFunctionBegin;
  info = TaoNew(TAO_NelderMead,&nm); CHKERRQ(info);
  tao->data = (void*)nm;

  info = PetscLogObjectMemory(tao,sizeof(TAO_NelderMead)); CHKERRQ(info);

  info = TaoSetTaoSolveRoutine(tao,TaoSolve_NelderMead,(void*)nm); CHKERRQ(info);
  info = TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_NelderMead, 
				    TaoSetDown_NelderMead); CHKERRQ(info);
  info = TaoSetTaoOptionsRoutine(tao,TaoSetOptions_NelderMead); CHKERRQ(info);
  info = TaoSetTaoViewRoutine(tao, TaoView_NelderMead); CHKERRQ(info);
  info = TaoSetMaximumIterates(tao,2000); CHKERRQ(info);
  info = TaoSetMaximumFunctionEvaluations(tao,4000); CHKERRQ(info);
  info = TaoSetTolerances(tao,1e-8,1e-8,0,0); CHKERRQ(info);


  nm->simplex = 0;
  nm->lamda = 1;

  nm->mu_ic = -0.5;
  nm->mu_oc = 0.5;
  nm->mu_r = 1.0;
  nm->mu_e = 2.0;

  TaoFunctionReturn(0);
}
EXTERN_C_END
    

/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "NelderMeadSort"
int NelderMeadSort(TAO_NelderMead *nm) {
  double *values = nm->f_values;
  int *indices = nm->indices;
  int dim = nm->N+1;

  int i,j,index;
  double val;
  TaoFunctionBegin;  
  for (i=1;i<dim;i++) {
    index = indices[i];
    val = values[index];
    for (j=i-1; j>=0 && values[indices[j]] > val; j--) {
      indices[j+1] = indices[j];
    }
    indices[j+1] = index;
  }  
  TaoFunctionReturn(0);
}


/*------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "NelderMeadReplace"
int NelderMeadReplace(TAO_NelderMead *nm, int index, TaoVec *Xmu, double f)
{
  int info;
  TaoFunctionBegin;
  // Add new vector's fraction of average
  info = nm->Xbar->Axpy(nm->oneOverN,Xmu); CHKERRQ(info);
  info = nm->simplex[index]->CopyFrom(Xmu); CHKERRQ(info);
  nm->f_values[index] = f;

  info = NelderMeadSort(nm);
  // Subtract last vector from average
  info = nm->Xbar->Axpy(-nm->oneOverN,
			nm->simplex[nm->indices[nm->N]]); CHKERRQ(info);
  TaoFunctionReturn(0);
}
