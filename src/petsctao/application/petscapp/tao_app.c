#include "tao_general.h"
#include "tao_app_impl.h"       /*I  "tao.h"  I*/

int Tao_ObjectiveEval=0, Tao_GradientEval=0, Tao_HessianEval=0,Tao_JacobianEval=0,Tao_FunctionEval=0;
int TAO_APP_COOKIE=0;

extern int TaoAppAddFiniteDifferences(TAO_APPLICATION);
extern int TaoAppDestroy(TAO_APPLICATION);

#undef __FUNCT__  
#define __FUNCT__ "TaoApplicationCreate"
/*@C
  TaoApplicationCreate - Creates a TaoApplication that
uses PETSc data structures.   The vectors used for gradient
and other information can be a PETSc Vec.  The routines
for function evaluation, and derivative information can
also used PETSc arguments.

   Input Parameters:
.  comm - an MPI communiicator

   Output Parameters:
.  newapp - the TaoApplication structure

.seealso TaoAppSetObjectiveAndGradientRoutine(), TaoSolveApplication(), TaoAppDestroy()

   Level: beginner

.keywords: Application
@*/
int TaoApplicationCreate(MPI_Comm comm, TAO_APPLICATION* newapp){
  int info;
  TAO_APPLICATION taoapp;

  PetscFunctionBegin;

  if (TAO_APP_COOKIE==0){
    info=PetscCookieRegister("TAO Application",&TAO_APP_COOKIE); CHKERRQ(info);
  }

  info = PetscHeaderCreate(taoapp,_p_TAOAPPLICATION,int,TAO_APP_COOKIE,-1,"TAO APP",comm,TaoAppDestroy,0); CHKERRQ(info);

  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,0);

  taoapp->computeumfunction=0; taoapp->computegradient=0;
  taoapp->computefunctiongradient=0; taoapp->computehessian=0;
  taoapp->computejacobian=0; taoapp->computevfunc=0;
  taoapp->computebounds=0; taoapp->boundctx=0;
  taoapp->grtol=0;
  taoapp->usrfctx=0; taoapp->usrfgctx=0; taoapp->usrgctx=0; taoapp->usrhctx=0;
  taoapp->V=0; 
  taoapp->G=0; taoapp->H=0; taoapp->HP=0; 
  taoapp->R=0; taoapp->J=0; taoapp->JP=0; 

  taoapp->numbermonitors     = 0;
  taoapp->numberdestroyers   = 0;
  taoapp->nAddOn             = 0;

  taoapp->numberoptioncheckers=0;

  if (Tao_ObjectiveEval==0){
    info = PetscLogEventRegister("TaoAppObjective",TAO_APP_COOKIE,&Tao_ObjectiveEval); CHKERRQ(info);
    info = PetscLogEventRegister("TaoAppGradient",TAO_APP_COOKIE,&Tao_GradientEval); CHKERRQ(info);
    info = PetscLogEventRegister("TaoAppHessian",TAO_APP_COOKIE,&Tao_HessianEval); CHKERRQ(info);
    info = PetscLogEventRegister("TaoAppJacobian",TAO_APP_COOKIE,&Tao_JacobianEval); CHKERRQ(info);
    info = PetscLogEventRegister("TaoAppFunction",TAO_APP_COOKIE,&Tao_FunctionEval); CHKERRQ(info);
  }
  info = TaoAppAddFiniteDifferences(taoapp); CHKERRQ(info);
  info = KSPCreate(comm,&taoapp->ksp); CHKERRQ(info);
  info = KSPSetFromOptions(taoapp->ksp); CHKERRQ(info);
  *newapp=taoapp;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppDestroy"
/*@
  TaoAppDestroy - Destroy the PETSc application
and all of the vectors and matrices associated with it.

   Input Parameters:
.  taoapp - the TaoApplication structure

.seealso TaoApplicationCreate(), TaoDestroy()

   Level: beginner

.keywords: Application, Destroy
@*/
int TaoAppDestroy(TAO_APPLICATION taoapp){
  int i,info;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (--((PetscObject)taoapp)->refct > 0) PetscFunctionReturn(0);
  if (taoapp->V){  info = VecDestroy(taoapp->V);  CHKERRQ(info); }
  if (taoapp->G){  info = VecDestroy(taoapp->G);  CHKERRQ(info); }
  if (taoapp->H){  info = MatDestroy(taoapp->H);  CHKERRQ(info); }
  if (taoapp->HP){ info = MatDestroy(taoapp->HP); CHKERRQ(info); }
  if (taoapp->R){  info = VecDestroy(taoapp->R);  CHKERRQ(info); }
  if (taoapp->J){  info = MatDestroy(taoapp->J);  CHKERRQ(info); }
  if (taoapp->JP){ info = MatDestroy(taoapp->JP); CHKERRQ(info); }
  
  if (taoapp->ksp) {
    info = KSPDestroy(taoapp->ksp); CHKERRQ(info);
    taoapp->ksp=0;
  }
  for (i=0; i< taoapp->numberdestroyers; i++){
    info = (*taoapp->userdestroy[i])(taoapp->userctxdestroy[i]); CHKERRQ(info);
  }
  
  //  delete taoapp->taoappl;
  
  PetscHeaderDestroy(taoapp); 
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetInitialSolutionVec"
/*@
   TaoAppSetInitialSolutionVec - Sets the vector representing the variables
   and an initial guess.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
-  xx - variable vector that stores the solution

   Level: beginner

   Note: 
   This vector will be used by the solver, so do not use it
   for other purposes.  The user should destroy this vector after
   solving the application.

   Note:  
   If the user is unaware of a decent initial solution,
   the vector should be set to zero.

   Note:  
   The TAO solvers will not use the contents of this 
   Vec until the TaoSolve() is called.  Therefore the user
   may compute an initial solution in this vector after this
   routine -- but before TaoSolve().

.seealso:  TaoAppGetSolutionVec(), TaoAppSetObjectiveRoutine()
@*/
int TaoAppSetInitialSolutionVec(TAO_APPLICATION taoapp, Vec xx){
  int info;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (xx){
    PetscValidHeaderSpecific(xx,VEC_COOKIE,2);
    PetscObjectReference((PetscObject)xx);
 }
  if (taoapp->V){
    info=VecDestroy(taoapp->V);CHKERRQ(info); 
  }  
  taoapp->V=xx;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetDefaultSolutionVec"
/*@
   TaoAppSetDefaultSolutionVec - Sets the vector representing the variables
   and an initial guess.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION context
-  xx - variable vector that stores the solution

   Level: beginner

   Note: 
   This vector will be used by the solver, so do not use it
   for other purposes.  The user should destroy this vector after
   solving the application.

.seealso:  TaoAppGetSolutionVec(), TaoAppSetObjectiveRoutine(), TaoAppSetInitialSolutionVec()
@*/
int TaoAppSetDefaultSolutionVec(TAO_APPLICATION taoapp, Vec xx){
  int info;
  PetscScalar zero=0.0;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (xx){
    PetscValidHeaderSpecific(xx,VEC_COOKIE,2);
    PetscObjectReference((PetscObject)xx);
 }
  if (taoapp->V){
    info=VecDestroy(taoapp->V);CHKERRQ(info); 
  }  
  taoapp->V=xx;
  info = VecSet(xx, zero); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetSolutionVec"
/*@
  TaoAppGetSolutionVec - Get the vector with the
  solution in the current application.

   Input Parameters:
.  taoapp - the application

   Output Parameter:
.  X - the solution vector

   Note: 
   This vector should not be destroyed.

   Level: intermediate


.keywords: Application, variables
@*/
int TaoAppGetSolutionVec(TAO_APPLICATION taoapp, Vec *X){
  PetscFunctionBegin; 
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (X){
    *X=taoapp->V;
  }
  PetscFunctionReturn(0);
}


/* ------------ Routines to set performance monitoring options ----------- */

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetMonitor"
/*@C
   TaoAppSetMonitor - Sets an ADDITIONAL function that is to be used at every
   iteration of the solver to display the iteration's 
   progress.   

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION solver context
.  mymonitor - monitoring routine
-  mctx - [optional] user-defined context for private data for the 
          monitor routine (may be TAO_NULL)

   Calling sequence of mymonitor:
$     int mymonitor(TAO_APPLICATION taoapp,void *mctx)

+    taoapp - the TAO_APPLICATION solver context
-    mctx - [optional] monitoring context


   Note: 
   Several different monitoring routines may be set by calling
   TaoAppSetMonitor() multiple times; all will be called in the 
   order in which they were set.

   Level: intermediate

.keywords: options, monitor, View

.seealso: TaoSetMonitor(), TaoAppSetDestroyRoutine()
@*/
int TaoAppSetMonitor(TAO_APPLICATION taoapp,int (*mymonitor)(TAO_APPLICATION,void*),void *mctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (mymonitor){
    if (taoapp->numbermonitors >= MAX_TAO_MONITORS) {
      SETERRQ(1,"Too many monitors set");
    }    
    taoapp->monitor[taoapp->numbermonitors]           = mymonitor;
    taoapp->monitorcontext[taoapp->numbermonitors++]  = (void*)mctx;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppMonitor"
/*@
   TaoAppMonitor - Apply the monitor functions for a TAO_APPLICATION object.

   Collective on TAO_APPLICATION

   Input Parameters:
.  taoapp - the TAO_APPLICATION structure

   Level: developer

.keywords: options, monitor, View

.seealso: TaoAppSetMonitor()
@*/
int TaoAppMonitor(TAO_APPLICATION taoapp){
  int i,info;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  for ( i=0; i<taoapp->numbermonitors; i++ ) {
    info = (*taoapp->monitor[i])(taoapp,taoapp->monitorcontext[i]);CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}


/* ------------ Routines to called when destroying this application ----------- */
#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetDestroyRoutine"
/*@C
   TaoAppSetDestroyRoutine - Sets an ADDITIONAL function that will be called when
   this application is destroyed.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION solver context
.  destroy - function pointer
-  ctx - [optional] user-defined context for private data for the 
          destroy routine (may be TAO_NULL)

   Calling sequence of destroy:
$     int mydestroy(void *ctx)

.    ctx - [optional] destroy context


   Level: intermediate

   Note:  
   This routine is often used to destroy structures used by monitors and 
   function evaluations.  This routine may also be used to shut down other packages
   such as ADIC.

.keywords: destroy

.seealso: TaoAppSetMonitor(), TaoAppSetHessianRoutine(), TaoAppDestroy()
@*/
int TaoAppSetDestroyRoutine(TAO_APPLICATION taoapp,int (*destroy)(void*),void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (destroy){
    if (taoapp->numberdestroyers >= MAX_TAO_USER_DESTROY) {
      SETERRQ(1,"TAO ERRROR: Too many TAO APPLICATION destroy routines set");
    }
    
    taoapp->userdestroy[taoapp->numberdestroyers]           = destroy;
    taoapp->userctxdestroy[taoapp->numberdestroyers++]      = ctx;
  }
  PetscFunctionReturn(0);
}

/* ------------ Routines to extend TaoApp ----------- */

#undef __FUNCT__  
#define __FUNCT__ "TaoAppAddObject"
/*@
   TaoAppAddObject -  add an object from to the Tao Application.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION solver context
.  key - string used to ID this object
-  ctx - user-defined context for private data

   Output Paramter:
.  id - ignored
   

   Note: 
   This routine can be used to extend the functionality of this object

   Level: advanced

.keywords: extensions

.seealso: TaoAppQueryForObject()
@*/
int TaoAppAddObject(TAO_APPLICATION taoapp, char *key, void *ctx, TaoInt *id)
{
  int info;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (ctx){
    if (taoapp->nAddOn >= MAX_TAOAPP_ADDON) {
      SETERRQ(1,"Too many TaoObject added on");
    }
    taoapp->TaoAppCtx[taoapp->nAddOn].ctx = ctx;
    info=PetscStrncpy(taoapp->TaoAppCtx[taoapp->nAddOn].key , key, MAX_TAO_KEY_LENGTH); CHKERRQ(info);
    taoapp->TaoAppCtx[taoapp->nAddOn].id = taoapp->nAddOn;
    *id=taoapp->TaoAppCtx[taoapp->nAddOn].id;
    taoapp->nAddOn++;
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppQueryForObject"
/*@
   TaoAppQueryForObject -  query the TAO Application for an object

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION solver context
.  key - string used to ID this object
-  ctx - user-defined context for private data

   Output Parameter:
   

   Note: 
   This routine can be used to extend the functionality of this object

   Level: advanced

.keywords: extensions

.seealso: TaoAppAddObject()
@*/
int TaoAppQueryForObject(TAO_APPLICATION taoapp, char *key, void **ctx)
{
  int i,n,info;
  PetscTruth flag=PETSC_FALSE;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  n=taoapp->nAddOn;
  *ctx=0;
  for (i=0;i<n;i++){
    info=PetscStrncmp(taoapp->TaoAppCtx[i].key,key,MAX_TAOAPP_ADDON,&flag); CHKERRQ(info);
    if (flag==PETSC_TRUE){
      *ctx=taoapp->TaoAppCtx[i].ctx;
      break;
    }
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "TaoAppQueryRemoveObject"
/*@
   TaoAppQueryRemoveObject -  add an object from to the Tao Application.

   Collective on TAO_APPLICATION

   Input Parameters:
+  taoapp - the TAO_APPLICATION solver context
-  key - string used to ID this object

   Note: 
   This routine can be used to extend the functionality of this object

   Level: advanced

.keywords: extensions

.seealso: TaoAppQueryForObject()
@*/
int TaoAppQueryRemoveObject(TAO_APPLICATION taoapp, char *key)
{
  int i,n,info;
  PetscTruth flag;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  n=taoapp->nAddOn;
  for (i=0;i<n;i++){
    info=PetscStrncmp(taoapp->TaoAppCtx[i].key,key,MAX_TAO_KEY_LENGTH,&flag); CHKERRQ(info);
    if (flag==PETSC_TRUE){
      if (n>0){
	taoapp->TaoAppCtx[i].ctx=taoapp->TaoAppCtx[n-1].ctx;
	taoapp->TaoAppCtx[i].id=taoapp->TaoAppCtx[n-1].id;
	info=PetscStrncpy(taoapp->TaoAppCtx[i].key,taoapp->TaoAppCtx[n-1].key,MAX_TAO_KEY_LENGTH); 
	CHKERRQ(info);
      }
      taoapp->nAddOn--;
      break;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetOptionsRoutine"
/*@C
   TaoAppSetOptionsRoutine - Sets an ADDITIONAL function that is to be called
   during TaoSetFromOptions().

   Collective on TAO_APPLICATION

   Input Parameters:
+  tao - the TAO_APPLICATION solver context
-  options -  routine the checks options

   Calling sequence of options:
$     int myoptions(TAO_APPLICATION taoapp)

.    taoapp - the TAO_APPLICATION solver context


   Note: 
   Several different options routines may be set by calling
   this routine.

   Level: advanced

.keywords: options, monitor, View

.seealso: TaoAppSetFromOptions()
@*/
int TaoAppSetOptionsRoutine(TAO_APPLICATION taoapp,int (*options)(TAO_APPLICATION) )
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (options){
    if (taoapp->numberoptioncheckers >= MAX_TAO_MONITORS) {
      SETERRQ(1,"Too many options checkers set");
    }    
    taoapp->checkoptions[taoapp->numberoptioncheckers++] = options;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppSetFromOptions"
/*@
  TaoAppSetFromOptions - Sets various TAO parameters from user options

   Collective on TAO_APPLICATION

   Input Parameters:
.  taoapp - the TAO Application

   Level: beginner

.keywords:  options

.seealso: TaoSetFromOptions()

@*/
int TaoAppSetFromOptions(TAO_APPLICATION taoapp){
  int i,info;
  const char *prefix=0;
  MPI_Comm comm;
  PetscTruth flg1;
  PetscFunctionBegin;

  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  info = PetscObjectGetOptionsPrefix((PetscObject)taoapp,&prefix); CHKERRQ(info);
  info = PetscObjectGetComm((PetscObject)taoapp,&comm);CHKERRQ(info);
  info = PetscOptionsBegin(comm,prefix,"TAO PETSC APPLICATIONS ","solver");CHKERRQ(info);
  info = PetscOptionsReal("-taoapp_grtol","the relative tolerance","TaoAppSetRelativeTolerance",
			  taoapp->grtol,&taoapp->grtol,&flg1);CHKERRQ(info);
  for (i=0;i<taoapp->numberoptioncheckers;i++){
    info = (*taoapp->checkoptions[i])(taoapp);CHKERRQ(info);
  }
  info = KSPSetFromOptions(taoapp->ksp);CHKERRQ(info);
  info = PetscOptionsEnd(); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "TaoAppGetKSP"
int TaoAppGetKSP(TAO_APPLICATION taoapp, KSP *ksp){
/*@C
  TaoAppGetKSP - Gets the linear solver used by the optimization application

   Collective on TAO_APPLICATION

   Input Parameters:
.  taoapp - the TAO Application

   Level: intermediate

.keywords:  options

.seealso: TaoSetFromOptions()

@*/
  PetscFunctionBegin; 
  PetscValidHeaderSpecific(taoapp,TAO_APP_COOKIE,1);
  if (ksp){
    *ksp=taoapp->ksp;
  }
  PetscFunctionReturn(0);
}
