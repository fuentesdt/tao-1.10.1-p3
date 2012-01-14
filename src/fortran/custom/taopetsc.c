#include "private/fortranimpl.h"
#include "petscksp.h"
#include "tao.h"


#ifdef PETSC_HAVE_FORTRAN_CAPS

#define taopetscapplicationcreate_     TAOPETSCAPPLICATIONCREATE
#define taoapplicationcreate_          TAOAPPLICATIONCREATE
#define taoapplicationdestroy_          TAOAPPLICATIONDESTROY
#define taosetapplication_          TAOSETAPPLICATION

#define taoappsetobjectiveandgradientroutine_    TAOAPPSETOBJECTIVEANDGRADIENTROUTINE
#define taoappsetobjectiveandgradientro_         TAOAPPSETOBJECTIVEANDGRADIENTRO
#define taoappsetobjectiveroutine_               TAOAPPSETOBJECTIVEROUTINE
#define taoappsetgradientroutine_                TAOAPPSETGRADIENTROUTINE
#define taoappsethessianroutine_                 TAOAPPSETHESSIANROUTINE
#define taoappsetvariableboundsroutine_          TAOAPPSETVARIABLEBOUNDSROUTINE
#define taoappsetjacobianroutine_                TAOAPPSETJACOBIANROUTINE
#define taoappsetconstraintroutine_              TAOAPPSETCONSTRAINTROUTINE
#define taoappsetmonitor_                        TAOAPPSETMONITOR
#define taoappgetksp_                            TAOAPPGETKSP

/* Grid application */
#define daappsethessianroutine_                         DAAPPSETHESSIANROUTINE
#define daappsetobjectiveandgradientroutine_            DAAPPSETOBJECTIVEANDGRADIENTROUTINE
#define daappsetobjectiveandgradientrou_                DAAPPSETOBJECTIVEANDGRADIENTROU
#define daappsetgradientroutine_                        DAAPPSETGRADIENTROUTINE
#define daappsetobjectiveroutine_                       DAAPPSETOBJECTIVEROUTINE
#define daappsetvariableboundsroutine_                  DAAPPSETVARIABLEBOUNDSROUTINE
#define daappsetconstraintroutine_                      DAAPPSETCONSTRAINTROUTINE
#define daappsetjacobianroutine_                        DAAPPSETJACOBIANROUTINE
#define daappsetbeforemonitor_                          DAAPPSETBEFOREMONITOR
#define daappsetaftermonitor_                           DAAPPSETAFTERMONITOR
#define daappsetelementobjectiveandgradientroutine_     DAAPPSETELEMENTOBJECTIVEANDGRADIENTROUTINE
#define daappsetelementobjectiveandgrad_                DAAPPSETELEMENTOBJECTIVEANDGRAD
#define daappsetelementhessianroutine_                  DAAPPSETELEMENTHESSIANROUTINE


#define taodefaultcomputehessian_   TAODEFAULTCOMPUTEHESSIAN
#define taodefaultcomputehessiancolor_ TAODEFAULTCOMPUTEHESSIANCOLOR

#define taogetksp_                      TAOGETKSP

#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)

#define taopetscapplicationcreate_     taopetscapplicationcreate 
#define taoapplicationcreate_                taoapplicationcreate 
#define taoapplicationdestroy_                taoapplicationdestroy 
#define taosetapplication_                taosetapplication 


#define taoappsetobjectiveandgradientroutine_    taoappsetobjectiveandgradientroutine
#define taoappsetobjectiveandgradientro_         taoappsetobjectiveandgradientro
#define taoappsetobjectiveroutine_               taoappsetobjectiveroutine
#define taoappsetgradientroutine_                taoappsetgradientroutine
#define taoappsethessianroutine_                 taoappsethessianroutine
#define taoappsetvariableboundsroutine_          taoappsetvariableboundsroutine
#define taoappsetjacobianroutine_                taoappsetjacobianroutine
#define taoappsetconstraintroutine_              taoappsetconstraintroutine
#define taoappsetmonitor_                        taoappsetmonitor
#define taoappgetksp_                            taoappgetksp


/* Grid application */
#define daappsethessianroutine_                         daappsethessianroutine
#define daappsetobjectiveandgradientroutine_            daappsetobjectiveandgradientroutine
#define daappsetobjectiveandgradientrou_                daappsetobjectiveandgradientrou
#define daappsetgradientroutine_                        daappsetgradientroutine
#define daappsetobjectiveroutine_                       daappsetobjectiveroutine
#define daappsetvariableboundsroutine_                  daappsetvariableboundsroutine
#define daappsetconstraintroutine_                      daappsetconstraintroutine
#define daappsetjacobianroutine_                        daappsetjacobianroutine
#define daappsetbeforemonitor_                          daappsetbeforemonitor
#define daappsetaftermonitor_                           daappsetaftermonitor
#define daappsetelementobjectiveandgradientroutine_     daappsetelementobjectiveandgradientroutine
#define daappsetelementobjectiveandgrad_                daappsetelementobjectiveandgrad
#define daappsetelementhessianroutine_                  daappsetelementhessianroutine

#define taodefaultcomputehessian_   taodefaultcomputehessian
#define taodefaultcomputehessiancolor_ taodefaultcomputehessiancolor

#define taogetksp_                      taogetksp

#endif

EXTERN_C_BEGIN

void PETSC_STDCALL taopetscapplicationcreate_(MPI_Comm *comm,TAO_APPLICATION  *outtao,int *ierr){

  *ierr = TaoPetscApplicationCreate(MPI_Comm_f2c(*(MPI_Fint *)&*comm),outtao);

}

void PETSC_STDCALL taoapplicationcreate_(MPI_Comm *comm,TAO_APPLICATION  *outtao,int *ierr){

  *ierr = TaoApplicationCreate(MPI_Comm_f2c(*(MPI_Fint *)&*comm),outtao);

}

void PETSC_STDCALL taoapplicationdestroy_(TAO_APPLICATION  *outtao,int *ierr){

  *ierr = TaoApplicationDestroy(*outtao);

}

/* -------------- Setting call-back routines ---------------- */


EXTERN_C_END

static int ourtaominfunctiongradientroutine(TAO_APPLICATION taoapp,Vec x,double *f,Vec r,void *ctx)
{
  int ierr = 0;
  (*(void (PETSC_STDCALL *)(TAO_APPLICATION*,Vec*,double*,Vec*,void*,int*))(((PetscObject)taoapp)->fortran_func_pointers[1]))(&taoapp,&x,f,&r,ctx,&ierr);
  return 0;
}

EXTERN_C_BEGIN


void PETSC_STDCALL taoappsetobjectiveandgradientroutine_(TAO_APPLICATION  *taoapp,void (PETSC_STDCALL *func)(TAO_APPLICATION*,Vec*,double*,Vec*,void*,int*),void *ctx,int *ierr){
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  ((PetscObject)*taoapp)->fortran_func_pointers[1] = (PetscVoidFunction)func;
  *ierr = TaoAppSetObjectiveAndGradientRoutine(*taoapp,ourtaominfunctiongradientroutine,ctx);
}
void PETSC_STDCALL taoappsetobjectiveandgradientro_(TAO_APPLICATION  *taoapp,void (PETSC_STDCALL *func)(TAO_APPLICATION*,Vec*,double*,Vec*,void*,int*),void *ctx,int *ierr){
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  ((PetscObject)*taoapp)->fortran_func_pointers[1] = (PetscVoidFunction)func;
  *ierr = TaoAppSetObjectiveAndGradientRoutine(*taoapp,ourtaominfunctiongradientroutine,ctx);
}


EXTERN_C_END
static int ourtaominfunction(TAO_APPLICATION taoapp,Vec x,double* d,void *ctx)
{
  int ierr = 0;
  (*(void (PETSC_STDCALL *)(TAO_APPLICATION*,Vec*,double*,void*,int*))(((PetscObject)taoapp)->fortran_func_pointers[2]))(&taoapp,&x,d,ctx,&ierr);CHKERRQ(ierr);
  return 0;
}
EXTERN_C_BEGIN


void PETSC_STDCALL taoappsetobjectiveroutine_(TAO_APPLICATION  *taoapp,
	  void (PETSC_STDCALL *func)(TAO_APPLICATION*,Vec*,double*,void*,int*),void *ctx,int *ierr){
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  ((PetscObject)*taoapp)->fortran_func_pointers[2] = (PetscVoidFunction)func;
  *ierr = TaoAppSetObjectiveRoutine(*taoapp,ourtaominfunction,ctx);
}



EXTERN_C_END

static int ourtaogradientfunction(TAO_APPLICATION taoapp,Vec x,Vec d,void *ctx)
{
  int ierr = 0;
  (*(void (PETSC_STDCALL *)(TAO_APPLICATION*,Vec*,Vec*,void*,int*))(((PetscObject)taoapp)->fortran_func_pointers[3]))(&taoapp,&x,&d,ctx,&ierr);CHKERRQ(ierr);
  return 0;
}
EXTERN_C_BEGIN


void PETSC_STDCALL taoappsetgradientroutine_(TAO_APPLICATION *taoapp,void (PETSC_STDCALL *func)(TAO_APPLICATION*,Vec*,Vec*,void*,int*),void *ctx,int *ierr){
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  ((PetscObject)*taoapp)->fortran_func_pointers[3] = (PetscVoidFunction)func;
  *ierr = TaoAppSetGradientRoutine(*taoapp,ourtaogradientfunction,ctx);
}


/* ---------------------------------------------------------*/
/*
  taodefaultcomputehessian() and taodefaultcomputehessiancolor().
  These routines can be used directly from Fortran, but the following is done
  primarily so that a Fortran call to TaoAppSetPetscHessian() will properly handle the
  defaults being passed in.

  functions, hence no STDCALL
*/
void taodefaultcomputehessian_(TAO_APPLICATION *tao,Vec *x,Mat *m,Mat *p,MatStructure* type,
                               void *ctx,int *ierr)
{
  *ierr = TaoAppDefaultComputeHessian(*tao,*x,m,p,type,ctx);
}

void  taodefaultcomputehessiancolor_(TAO_APPLICATION *tao,Vec *x,Mat *m,Mat *p,
                                     MatStructure* type,void *ctx,int *ierr)
{
  *ierr = TaoAppDefaultComputeHessianColor(*tao,*x,m,p,type,*(MatFDColoring*)ctx);
}

EXTERN_C_END

static int ourtaohessian(TAO_APPLICATION taoapp,Vec x,Mat* m,Mat* p,MatStructure* type,void*ctx)
{
  int ierr = 0;
  (*(void (PETSC_STDCALL *)(TAO_APPLICATION*,Vec*,Mat*,Mat*,MatStructure*,void*,int*))(((PetscObject)taoapp)->fortran_func_pointers[4]))(&taoapp,&x,m,p,type,ctx,&ierr);CHKERRQ(ierr);
  return 0;
}

EXTERN_C_BEGIN


void PETSC_STDCALL taoappsethessianroutine_(TAO_APPLICATION  *taoapp,void (PETSC_STDCALL *func)(TAO_APPLICATION*,Vec*,Mat*,Mat*,
                                  MatStructure*,void*,int*),void *ctx,int *ierr)
{
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  if ((void(*)())func == (void(*)()) taodefaultcomputehessian_) {
    *ierr = TaoAppSetHessianRoutine(*taoapp,TaoAppDefaultComputeHessian,ctx);
  } else if ( (void(*)()) func == (void(*)()) taodefaultcomputehessiancolor_) {
    *ierr = TaoAppSetHessianRoutine(*taoapp,TaoAppDefaultComputeHessianColor,*(MatFDColoring*)ctx);
  } else {

  ((PetscObject)*taoapp)->fortran_func_pointers[4] = (PetscVoidFunction)func;
    *ierr = TaoAppSetHessianRoutine(*taoapp,ourtaohessian,ctx);
  }
}



EXTERN_C_END

static int ourtaovariableboundsroutine(TAO_APPLICATION taoapp, Vec x,Vec d,void *ctx)
{
  int ierr = 0;
  (*(void (PETSC_STDCALL *)(TAO_APPLICATION*,Vec*,Vec*,void*,int*))(((PetscObject)taoapp)->fortran_func_pointers[5]))(&taoapp,&x,&d,ctx,&ierr);CHKERRQ(ierr);
  return 0;
}

EXTERN_C_BEGIN


void PETSC_STDCALL taoappsetvariableboundsroutine_(TAO_APPLICATION *taoapp,void (PETSC_STDCALL *func)(TAO_APPLICATION*,Vec*,Vec*,void*,int*),void *ctx,int *ierr){
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  ((PetscObject)*taoapp)->fortran_func_pointers[5] = (PetscVoidFunction)func;
  *ierr = TaoAppSetVariableBoundsRoutine(*taoapp,ourtaovariableboundsroutine,ctx);
}


EXTERN_C_END
static int ourtaoconstraintsfunction(TAO_APPLICATION taoapp,Vec x,Vec d,void *ctx)
{
  int ierr = 0;
  (*(void (PETSC_STDCALL *)(TAO_APPLICATION*,Vec*,Vec*,void*,int*))(((PetscObject)taoapp)->fortran_func_pointers[6]))(&taoapp,&x,&d,ctx,&ierr);CHKERRQ(ierr);
  return 0;
}

static int ourtaoappsetmonitor(TAO_APPLICATION taoapp, void *ctx) {
  int ierr = 0;
  (*(void (PETSC_STDCALL *)(TAO_APPLICATION*,void*,int*))(((PetscObject)taoapp)->fortran_func_pointers[7]))(&taoapp,ctx,&ierr); CHKERRQ(ierr);
  return 0;
}
EXTERN_C_BEGIN


void PETSC_STDCALL taoappsetconstraintroutine_(TAO_APPLICATION  *taoapp,void (PETSC_STDCALL *func)(TAO_APPLICATION*,Vec*,Vec*,void*,int*),void *ctx,int *ierr){
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  ((PetscObject)*taoapp)->fortran_func_pointers[6] = (PetscVoidFunction)func;
  *ierr = TaoAppSetConstraintRoutine(*taoapp,ourtaoconstraintsfunction,ctx);
}

void PETSC_STDCALL taoappsetmonitor_(TAO_APPLICATION *taoapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*,void *, int*), void *ctx, int *ierr){
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  ((PetscObject)*taoapp)->fortran_func_pointers[7] = (PetscVoidFunction)func;
  *ierr = TaoAppSetMonitor(*taoapp,ourtaoappsetmonitor,ctx);
}


EXTERN_C_END

static int ourtaojacobianfunction(TAO_APPLICATION taoapp,Vec x,Mat *J,Mat *JP, MatStructure *flag,void *ctx)
{
  int ierr = 0;
  (*(void (PETSC_STDCALL *)(TAO_APPLICATION*,Vec*,Mat*,Mat*,MatStructure*,void*,int*))(((PetscObject)taoapp)->fortran_func_pointers[8]))(&taoapp,&x,J,JP,flag,ctx,&ierr);CHKERRQ(ierr);
  return 0;
}
EXTERN_C_BEGIN


void PETSC_STDCALL taoappsetjacobianroutine_(TAO_APPLICATION  *taoapp,void (PETSC_STDCALL *func)(TAO_APPLICATION*,Vec*,Mat*,Mat*,MatStructure*,void*,int*),void *ctx,int *ierr){
  CHKFORTRANNULLFUNCTION(func);
  PetscObjectAllocateFortranPointers(*taoapp,12);
  ((PetscObject)*taoapp)->fortran_func_pointers[8] = (PetscVoidFunction)func;
  *ierr = TaoAppSetJacobianRoutine(*taoapp,ourtaojacobianfunction,ctx);
}
void PETSC_STDCALL taoappgetksp_(TAO_APPLICATION *taoapp,KSP *tksp, int *ierr ){
  *ierr = TaoAppGetKSP(*taoapp,tksp);
}


void PETSC_STDCALL taogetksp_(TAO_SOLVER  *tao,KSP *tksp, int *ierr ){
  *ierr = TaoGetKSP(*tao,tksp);
}


EXTERN_C_END

#include "taodaapplication.h"

/* Grid Application routines */
EXTERN_C_BEGIN
static void (PETSC_STDCALL *grid1)(TAO_APPLICATION*, DA*, Vec*, Mat*, void*, int*);
static void (PETSC_STDCALL *grid2)(TAO_APPLICATION*, DA*, Vec*, double*, Vec*, void*, int*);
static void (PETSC_STDCALL *grid3)(TAO_APPLICATION*, DA*, Vec*, Vec*, void*, int*);
static void (PETSC_STDCALL *grid4)(TAO_APPLICATION*, DA*, Vec*, double*, void*, int*);
static void (PETSC_STDCALL *grid5)(TAO_APPLICATION*, DA*, Vec*, Vec*, void*, int*);
static void (PETSC_STDCALL *grid6)(TAO_APPLICATION*, DA*, Vec*, Vec*, void*, int*);
static void (PETSC_STDCALL *grid7)(TAO_APPLICATION*, DA*, Vec*, Mat*, void*, int*);
static void (PETSC_STDCALL *grid8)(TAO_APPLICATION*, DA*, TaoInt*, void*, int*);
static void (PETSC_STDCALL *grid9)(TAO_APPLICATION*, DA*, TaoInt*, void*, int*);
/* Grid Application Element routines */
static void (PETSC_STDCALL *grid10)(TaoInt[2], double[4], double*, double[4], void*, int*);
static void (PETSC_STDCALL *grid11)(TaoInt[2], double[4], double[4][4], void*, int*);
EXTERN_C_END
static int ourdaappsethessianroutine(TAO_APPLICATION daapp, DA da, Vec x, Mat H, void *ctx) {
  int ierr = 0;
  (*grid1)(&daapp, &da, &x, &H, ctx, &ierr);
  return ierr;
}
static int ourdaappsetobjectiveandgradientroutine(TAO_APPLICATION daapp, DA da, Vec x, double *f, Vec G, void *ctx) {
  int ierr = 0;
  (*grid2)(&daapp, &da, &x, f, &G, ctx, &ierr);
  return ierr;
}
static int ourdaappsetgradientroutine(TAO_APPLICATION daapp, DA da, Vec x, Vec g, void *ctx) {
  int ierr = 0;
  (*grid3)(&daapp, &da, &x, &g, ctx, &ierr);
  return ierr;
}
static int ourdaappsetobjectiveroutine(TAO_APPLICATION daapp, DA da, Vec x, double *f, void *ctx) {
  int ierr = 0;
  (*grid4)(&daapp, &da, &x, f, ctx, &ierr);
  return ierr;
}
static int ourdaappsetvariableboundsroutine(TAO_APPLICATION daapp, DA da, Vec L, Vec U, void *ctx) {
  int ierr = 0;
  (*grid5)(&daapp, &da, &L, &U, ctx, &ierr);
  return ierr;
}
static int ourdaappsetconstraintroutine(TAO_APPLICATION daapp, DA da, Vec X, Vec R, void *ctx) {
  int ierr = 0;
  (*grid6)(&daapp, &da, &X, &R, ctx, &ierr);
  return ierr;
}
static int ourdaappsetjacobianroutine(TAO_APPLICATION daapp, DA da, Vec X, Mat J, void *ctx) {
  int ierr = 0;
  (*grid7)(&daapp, &da, &X, &J, ctx, &ierr);
  return ierr;
}
static int ourdaappsetbeforemonitor(TAO_APPLICATION daapp, DA da, TaoInt levels, void *ctx) {
  int ierr = 0;
  (*grid8)(&daapp, &da, &levels, ctx, &ierr);
  return ierr;
}
static int ourdaappsetaftermonitor(TAO_APPLICATION daapp, DA da, TaoInt levels, void *ctx) {
  int ierr = 0;
  (*grid9)(&daapp, &da, &levels, ctx, &ierr);
  return ierr;
}
static int ourdaappsetelementobjectiveandgradientroutine(PetscInt coord[2], double x[4], double *f, double g[4], void *ctx) {
  int ierr = 0;
  (*grid10)(coord, x, f, g, ctx, &ierr);
  return ierr;
}
static int ourdaappsetelementhessianroutine(PetscInt coord[2], double x[4], double H[4][4], void *ctx) {
  int ierr = 0;
  (*grid11)(coord, x, H, ctx, &ierr);
  return ierr;
}
EXTERN_C_BEGIN
void PETSC_STDCALL daappsethessianroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*,DA*,Vec*, Mat*, void*, int*),void *ctx, int *ierr) {
  grid1 = func;
  *ierr = DAAppSetHessianRoutine(*daapp, ourdaappsethessianroutine, ctx);
}
void PETSC_STDCALL daappsetobjectiveandgradientroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, Vec*, double*, Vec*, void*, int*), void *ctx, int *ierr) {
  grid2 = func;
  *ierr = DAAppSetObjectiveAndGradientRoutine(*daapp, ourdaappsetobjectiveandgradientroutine, ctx);
}
void PETSC_STDCALL daappsetobjectiveandgradientrou_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, Vec*, double*, Vec*, void*, int*), void *ctx, int *ierr) {
  grid2 = func;
  *ierr = DAAppSetObjectiveAndGradientRoutine(*daapp, ourdaappsetobjectiveandgradientroutine, ctx);
}

void PETSC_STDCALL daappsetgradientroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, Vec*, Vec*, void*, int*), void *ctx, int *ierr) {
  grid3 = func;
  *ierr = DAAppSetGradientRoutine(*daapp, ourdaappsetgradientroutine, ctx);
}
void PETSC_STDCALL daappsetobjectiveroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, Vec*, double*, void*, int*), void *ctx, int *ierr) {
  grid4 = func;
  *ierr = DAAppSetObjectiveRoutine(*daapp, ourdaappsetobjectiveroutine, ctx);
}
void PETSC_STDCALL daappsetvariableboundsroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, Vec*, Vec*, void*, int*), void *ctx, int *ierr) {
  grid5 = func;
  *ierr = DAAppSetVariableBoundsRoutine(*daapp, ourdaappsetvariableboundsroutine, ctx);
}
void PETSC_STDCALL daappsetconstraintroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, Vec*, Vec*, void*, int*), void *ctx, int *ierr) {
  grid6 = func;
  *ierr = DAAppSetConstraintRoutine(*daapp, ourdaappsetconstraintroutine, ctx);
}
void PETSC_STDCALL daappsetjacobianroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, Vec*, Mat*, void*, int*), void *ctx, int *ierr) {
  grid7 = func;
  *ierr = DAAppSetJacobianRoutine(*daapp, ourdaappsetjacobianroutine, ctx);
}
void PETSC_STDCALL daappsetbeforemonitor_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, PetscInt*, void*, int*), void *ctx, int *ierr) {
  grid8 = func;
  *ierr = DAAppSetBeforeMonitor(*daapp, ourdaappsetbeforemonitor, ctx);
}
void PETSC_STDCALL daappsetaftermonitor_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(TAO_APPLICATION*, DA*, PetscInt*, void*, int*), void *ctx, int *ierr) {
  grid9 = func;
  *ierr = DAAppSetAfterMonitor(*daapp, ourdaappsetaftermonitor, ctx);
}
void PETSC_STDCALL daappsetelementobjectiveandgradientroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(PetscInt[2], double[4], double*, double[4], void*, int*), PetscInt *flops, void *ctx, int *ierr) {
  grid10 = func;
  *ierr = DAAppSetElementObjectiveAndGradientRoutine(*daapp, ourdaappsetelementobjectiveandgradientroutine, *flops, ctx);
}
void PETSC_STDCALL daappsetelementobjectiveandgrad_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(PetscInt[2], double[4], double*, double[4], void*, int*), PetscInt *flops, void *ctx, int *ierr) {
  grid10 = func;
  *ierr = DAAppSetElementObjectiveAndGradientRoutine(*daapp, ourdaappsetelementobjectiveandgradientroutine, *flops, ctx);
}
void PETSC_STDCALL daappsetelementhessianroutine_(TAO_APPLICATION *daapp, void (PETSC_STDCALL *func)(PetscInt[2], double[4], double[4][4], void*, int*), PetscInt *flops, void *ctx, int *ierr) {
  grid11 = func;
  *ierr = DAAppSetElementHessianRoutine(*daapp, ourdaappsetelementhessianroutine, *flops, ctx);
}
EXTERN_C_END
