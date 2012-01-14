#include "tao_solver.h"
#include "src/tao_impl.h"
#include "apps_solver.h"
#include "apps.H"
#include "gci.H"
#include "args.H"
#include "fevalwkr.H"
#include "cachewkr.H"
#include "basepoint.H"
#include "taofevalwkr.h"
#include "taofevalmgr.h"






#undef __FUNCT__
#define __FUNCT__ "TaoSetUp_APPS"
int TaoSetUp_APPS(TAO_SOLVER tao, void *solver)
{
  int info;
  TAO_APPS *appsPtr = (TAO_APPS*)solver;
  TaoVec *xx;
  int size;

  TaoFunctionBegin;

  /* Need to make a gci communicator (duplicate the tao communicator) */
  // NOTICE: I had to make GCI::APPS_COMM public

  // info = MPI_Comm_dup(tao->comm,&GCI::APPS_COMM);


  
  info = TaoGetSolution(tao, &xx); CHKERRQ(info);
  info = xx->GetDimension(&appsPtr->ndim);

  // Set up static variables in the TaoFevalMgr class
  info = TaoFevalMgr::setTao(tao,appsPtr); CHKERRQ(info);


  info = TaoGetSolution(tao,&xx); CHKERRQ(info);
  info = xx->Clone(&appsPtr->xl); CHKERRQ(info);
  info = xx->Clone(&appsPtr->xu); CHKERRQ(info);
  info = TaoSetVariableBounds(tao,appsPtr->xl,appsPtr->xu);CHKERRQ(info);


  info = MPI_Comm_size(tao->comm,&size); CHKERRQ(info);
  // Check for enough processes
  if (size < 2*appsPtr->ndim + 2)
  {
    info = 1;  // loop reached
    SETERRQ(1,"Number of Processes must be at least 2*ndim + 2");
  }


  TaoFunctionReturn(info); 
}

//================================
#undef __FUNCT__
#define __FUNCT__ "TaoDestroy_APPS"
static int TaoDestroy_APPS(TAO_SOLVER tao, void *solver)
{

  TAO_APPS *appsPtr = (TAO_APPS*)solver;
  int i,info;

  TaoFunctionBegin;




  //free argv (These were used to pass options to APPS::main() 
  /*
  for (i=0;i<appsPtr->argc;i++)
  {
    info = TaoFree(appsPtr->argv[i]); CHKERRQ(info);
  }
  info = TaoFree(appsPtr->argv); CHKERRQ(info);
  */

  TaoFunctionReturn(0);
}


//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoSolve_APPS"
static int TaoSolve_APPS(TAO_SOLVER tao, void *solver)
{
  TaoVec *XX;
  double *x;
  int info;
  TAO_APPS *appsPtr = (TAO_APPS*)solver;
  FevalWkr *fevalwkrptr;
  int returnvalue;
  int taskid;
  int size;
  int i;
  int dim,nn;
  BasePoint *point;
  CacheWkr cachewkr;
  

  TaoFunctionBegin;

  /* Get options */
  
  
  /* Make sure that the local and global size of XX are the same */
  info = TaoGetSolution(tao, &XX); CHKERRQ(info);
  info = XX->GetDimension(&dim); CHKERRQ(info);
  info = XX->GetArray(&x,&nn); CHKERRQ(info);
  if (dim != nn)
  {
    SETERRQ(1,"TAO/APPS -- Local and global vector lengths must be equal");
  }
  info = XX->RestoreArray(&x,&nn); CHKERRQ(info);


  info = MPI_Comm_rank(tao->comm,&taskid); CHKERRQ(info);
  info = MPI_Comm_size(tao->comm,&size); CHKERRQ(info);

  if (taskid == 0)   // The main function is run on the rank 0 machine
  {

    // Call the apps solver, point is allocated and holds the minimizer
    point = APPS::main(appsPtr->argc, appsPtr->argv);

    if (point == NULL) {
      SETERRQ(1,"Failure to allocate point in APPS::main()");
    }
      

    info = XX->GetArray(&x,&nn); CHKERRQ(info);

    
    for (i=0;i<appsPtr->ndim;i++)
      x[i] = (*point)[i];

    info = XX->RestoreArray(&x,&nn); CHKERRQ(info);

    info = TaoComputeFunction(tao,XX,&appsPtr->fval); CHKERRQ(info);

    // Gather all the monitoring data into process 0
    info = gather(tao, appsPtr); CHKERRQ(info);


    // Deallocate point
    delete point;

    // These were intended to be delete from GCI::exit, so made them public and deleted here
    delete[] GCI::sendbuf;
    delete[] GCI::recvbuf;

    returnvalue = 0;
  }
  
  // The last machine runs the cachewkr (if selected)
  else if ((taskid == size - 1) && appsPtr->usecache)  
  {
    returnvalue = cachewkr.startWorking(); 
    gather(tao, appsPtr);
  }

  // Every other machine runs an fevalwkr
  else
  {
    fevalwkrptr = new TaoFevalWkr(tao,appsPtr);
    returnvalue = fevalwkrptr->startWorking();
    delete fevalwkrptr;
    info = gather(tao,appsPtr); CHKERRQ(info);
  }


  // Termination because step size is small
  tao->reason = TAO_CONVERGED_TRTOL;


  
  TaoFunctionReturn(returnvalue);

}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoSetOptions_APPS"
static int TaoSetOptions_APPS(TAO_SOLVER tao, void *solver)
{
  int info;
  TAO_APPS *appsPtr = (TAO_APPS*)solver;


  // Convert "-apps_option X" command line arguments into "--option=X" for use by APPS
  //  info = getoptions(appsPtr); CHKERRQ(info);

  TaoFunctionBegin;
  TaoFunctionReturn(0);
}

//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoView_APPS"
static int TaoView_APPS(TAO_SOLVER tao, void *solver)
{
  TaoFunctionBegin;
  TaoFunctionReturn(0);
}


//========================================
#undef __FUNCT__
#define __FUNCT__ "TaoCreate_APPS"
int TaoCreate_APPS(TAO_SOLVER tao)
{
  int info;
  TAO_APPS *appsPtr;

  TaoFunctionBegin;


  info = TaoNew(TAO_APPS, &appsPtr); CHKERRQ(info);

  info=TaoSetTaoSolveRoutine(tao,TaoSolve_APPS,(void*)appsPtr); CHKERRQ(info);
  info=TaoSetTaoSetUpDownRoutines(tao,TaoSetUp_APPS,TaoDestroy_APPS); CHKERRQ(info);
  info=TaoSetTaoOptionsRoutine(tao,TaoSetOptions_APPS); CHKERRQ(info);
  info=TaoSetTaoViewRoutine(tao,TaoView_APPS); CHKERRQ(info);

  appsPtr->usecache = 1;
  info = TaoSetMaximumIterates(tao,10000); CHKERRQ(info);

  GCI::APPS_COMM = tao->comm;

  TaoFunctionReturn(0);
}

//=======================================
#undef __FUNCT__
#define __FUNCT__ "getoptions"
int getoptions(TAO_APPS *appsPtr)
{
  // Haven't done break or reassign options
  // -o only works for PVM
  char options[16][10] = {"cache","alpha","pattern","search","tol","active","debug","profile",
			  "precision","inc","contract","step","min","max","nnls","eps"};
  char strings[35][50];  // Don't need this many, but just being safe
  char buffer[50]={0};
  char tempstring[16];
  int i;
  int info;
  PetscTruth flg,cacheisoff;
  int length;

  // Read in all options, store as strings.  Later, copy strings to argv


  appsPtr->argc=0;

  // don't really need this one, but need to start strings at argv[1]
  PetscStrcpy(strings[(appsPtr->argc)++],"executable"); 


  // It's possible to use PetscStrcpy and PetscStrcat instead of snprintf
    
  info = PetscOptionsGetString(0,"-apps_i",buffer, 49, &flg); CHKERRQ(info);
  if (flg)
  {
    snprintf(strings[(appsPtr->argc)++],50,"-i%-s",buffer);
  }

  info = PetscOptionsGetString(0,"-apps_o",buffer, 49, &flg); CHKERRQ(info);
  if (flg)
  {
    snprintf(strings[(appsPtr->argc)++],50,"-o%-s",buffer);
  }


  for (i=0;i<16;i++)
  {
    snprintf(tempstring,16,"-apps_%s",options[i]);
    info = PetscOptionsGetString(0,tempstring,buffer,49,&flg); CHKERRQ(info);
    if (flg)
    {
      snprintf(strings[(appsPtr->argc)++],50,"--%s=%-s",options[i],buffer);
      // The TaoFevalMgr needs to know if the cache is to be used
      // Assume it is to be used unless "-apps_cache false" is used
      if (i==0)
      {
	info = PetscStrncmp(buffer,"false",5,&cacheisoff); CHKERRQ(info);
	if (cacheisoff)
	  appsPtr->usecache = 0; 
      }
    }
  }


  info = TaoMalloc(sizeof(char *)*appsPtr->argc, &appsPtr->argv); CHKERRQ(info);
  
  for (i=0;i<appsPtr->argc; i++)
  {
    info = PetscStrlen(strings[i],&length); CHKERRQ(info);
    info = TaoMalloc(sizeof(char)*length+1, &appsPtr->argv[i]); CHKERRQ(info);
    info = PetscStrncpy(appsPtr->argv[i],strings[i],length); CHKERRQ(info);
  }


  TaoFunctionReturn(0);
}


//=======================================
#undef __FUNCT__
#define __FUNCT__ "gather"
int gather(TAO_SOLVER tao, TAO_APPS *appsPtr)
{
  double f;
  int local_nevaluations;
  int nevaluations=0;
  int info;
  int myid;
  double fnorm,cnorm,xdiff;
  TaoTerminateReason reason;


  TaoFunctionBegin;


  // The main process and the cachewkr process do not have anything to contribute
  myid = GCI::mytid();
  if ((myid == 0) || ((myid == GCI::getNumProcs()-1) && (appsPtr->usecache)))
  {
    local_nevaluations = 0;
  }

  else
  {
    info = TaoGetIterationData(tao, &local_nevaluations, &f, &fnorm, &cnorm, &xdiff, &reason); CHKERRQ(info);
  }

  // Combine the local nevaluations into nevaluations
  info = MPI_Reduce(&local_nevaluations, &nevaluations, 1, MPI_INT, MPI_SUM, 0, tao->comm); CHKERRQ(info);
  


  if (myid==0)
  {
    // Now report this information to the monitor
    info = TaoMonitor(tao, nevaluations, appsPtr->fval, 1, 0, 1, &reason); CHKERRQ(info);
  }

  TaoFunctionReturn(0);
  
  
}

int setargs(TAO_SOLVER tao, int argc, char *argv[])
{
  TAO_APPS *appsPtr = (TAO_APPS*)tao->data;
  appsPtr->argc = argc;
  appsPtr->argv = argv;
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "TaoRegisterAppspack"
int TaoRegisterAppspack(void)
{
  char path[] = "${TAO_DIR}/lib/${PETSC_ARCH}/libappspacktao.so";
  int info;
  
  TaoFunctionBegin;
  /* TODO:  why does path have to be zero? */
  info = TaoRegisterDynamic("tao_apps",0,"TaoCreate_APPS",TaoCreate_APPS); CHKERRQ(info);

  TaoFunctionReturn(0);
}





