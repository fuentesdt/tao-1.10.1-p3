/*$Id: minsurf3.c,v 1.7 2001/12/18 06:08:18 benson Exp $*/

/* Program usage: mpirun -np 1 rosenbrock1 [-help] [all TAO options] */

static  char help[] = "This example demonstrates use of the TAO package to \n\
solve an unconstrained minimization problem on a single processor.  We \n\
minimize the extended Rosenbrock function: \n\
   sum_{i=0}^{n/2-1} ( alpha*(x_{2i+1}-x_{2i}^2)^2 + (1-x_{2i})^2 ) \n";

/*T 
   Concepts: TAO - Solving an unconstrained minimization problem
   Routines: TaoInitialize(); TaoFinalize();
   Routines: TaoCreate(); TaoSetFunctionGradient(); 
   Routines: TaoSetHessian(); TaoSetInitialVector(); 
   Routines: TaoSolve(); TaoView(); TaoDestroy();
   Processors: 1
T*/ 

/* ---------------------------------------------------------------------- */

#include <iostream.h>
#include <math.h>
#include "Epetra_ESI.h"
#include "tao_solver.h"
#include "petscda.h"
#include "tao_esi.h"

/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines that evaluate the function,
   gradient, and hessian.
*/

class MinsurfExample: public esi::ESIApplication {

  
public:

  Vec        globalX,globalG,localX,orderX;
  DA         da;
  int        mx, my;
  int        bmx, bmy;
  double     bheight;
  double     *bottom, *top, *left, *right;             /* boundary values */

  ESIVCTR *xvec;
  esi::Operator<double,int>* H;
  epetra_esi::CrsMatrix<double,int>* petraH;


  MinsurfExample();
  ~MinsurfExample();

  int MSA_BoundaryConditions();
  int create_MESI_Operator(int numGlobal, int numLocal, esi::IndexSpace<int>*);

  virtual int getVariableVector(ESIVCTR **);
  virtual int getHessianMatrix(ESIMAT **HH);
  virtual int evaluateObjectiveAndGradient(ESIVCTR *, double *, ESIVCTR *);
  virtual int evaluateHessian(ESIVCTR *, ESIMAT *);
  virtual int monitor();

  int evaluateObjectiveFunction(ESIVCTR *, double *){return 1;}
  int evaluateGradient(ESIVCTR *xx, ESIVCTR *gg){return 1;}
  int evaluateBounds(esi::Vector<double,int> *,esi::Vector<double,int> *);

  virtual int initializeVariables(ESIVCTR *x){return 1;}


};

//typedef esi::Vector<double,int> ESIVCTR;

/* User-defined routines */


int cout_ESI_RowMatrix(esi::MatrixRowReadAccess<double,int>* esimat);
int cout_ESI_Vector(esi::Vector<double,int>* esivec);

int create_ESI_Solver(esi::Operator<double,int>* A,
                      esi::Solver<double,int>*& solver);

int create_ESI_IndexSpace(int numGlobal, int numLocal, Epetra_Comm& comm, 
			  esi::IndexSpace<int>*& esi_indexspace);

int create_ESI_Vector(esi::IndexSpace<int>* indexspace, esi::Vector<double,int>*& vec);

MinsurfExample::~MinsurfExample(){
  delete this->xvec;
  delete this->H;
  return;
}

MinsurfExample::MinsurfExample():esi::ESIApplication(){

  int i,j,info,low,high;
  int xs,ys,xm,ym,row;
  double *globalrow;
  PetscTruth flg;
  Epetra_Comm* comm = NULL;

  //  ErrorCode info;
  esi::Vector<double,int>* x,*g;

  int             Nx, Ny;              /* number of processors in x- and y- directions */
  int numLocal,numGlobal;
 int             size;


  this->mx=10;
  this->my=10;

#ifdef EPETRA_MPI
   comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
   comm = new Epetra_SerialComm();
#endif


  /*
    Create distributed array (DA) to manage parallel grid and vectors
  */
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  Nx = PETSC_DECIDE; Ny = PETSC_DECIDE;
  info = PetscOptionsGetInt(PETSC_NULL,"-mx",&this->mx,&flg);
  info = PetscOptionsGetInt(PETSC_NULL,"-my",&this->my,&flg);
  info = PetscOptionsGetInt(PETSC_NULL,"-Nx",&Nx,&flg);
  info = PetscOptionsGetInt(PETSC_NULL,"-Ny",&Ny,&flg);
  if (Nx*Ny != size && (Nx != PETSC_DECIDE || Ny != PETSC_DECIDE)){return;}
    //    SETERRQ(1,"Incompatible number of processors:  Nx * Ny != size");

  info = DACreate2d(MPI_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,this->mx,
                    this->my,Nx,Ny,1,1,PETSC_NULL,PETSC_NULL,&this->da);

  /*
     Extract global and local vectors from DA; The local vectors are
     used solely as work space for the evaluation of the function, 
     gradient, and Hessian.  Duplicate for remaining vectors that are 
     the same types.
  */
  info = DACreateGlobalVector(this->da,&this->globalX);
  info = VecDuplicate(this->globalX,&this->globalG);
  info = DACreateLocalVector(this->da,&this->localX);
  info = VecDuplicate(this->localX,&this->orderX);

  info = DAGetCorners(this->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);
  info = VecGetOwnershipRange(globalX,&low,&high);
  info = VecGetArray(globalX,&globalrow);

  for (j=0; j<ym; j++){
    for (i=0; i< xm; i++){    
      row=j*xm+i;
      globalrow[row]=low+row;
    }
  }

  info = VecRestoreArray(globalX,&globalrow);
  info = DAGlobalToLocalBegin(this->da,globalX,INSERT_VALUES,orderX);
  info = DAGlobalToLocalEnd(this->da,globalX,INSERT_VALUES,orderX);

  /* Calculate the boundary values */
  info = this->MSA_BoundaryConditions();

  info = VecGetLocalSize(this->globalX,&numLocal);
  info = VecGetSize(this->globalX,&numGlobal);

  /* Allocate vectors for the solution and gradient */
  esi::IndexSpace<int>* esiIndexSpace = NULL;
  info=create_ESI_IndexSpace(numGlobal, numLocal, *comm, esiIndexSpace); 
  //return error-code -1 if the allocation failed
  if (esiIndexSpace == NULL) return;

  MPI_Comm* mpicomm;
  esiIndexSpace->getRunTimeModel("MPI", (void*&)mpicomm);

  info = create_ESI_Vector(esiIndexSpace, x);

  if (x == NULL) return;

  /* 
     Create a matrix data structure to store the Hessian and set 
     the Hessian evalution routine.   
  */
  info = this->create_MESI_Operator(numGlobal, numLocal, esiIndexSpace);

  this->xvec=x;

  return;
}

int MinsurfExample::evaluateBounds(esi::Vector<double,int> *XL,esi::Vector<double,int> *XU){
  int i,j,row,info;
  int xs,ys,xm,ym;
  int mx=this->mx, my=this->my, bmy, bmx;
  double t1,t2,t3;
  double *xl;
  PetscTruth flg,cylinder;

  this->bmx=this->mx/2;
  this->bmy=this->my/2;
  this->bheight=0.3;
  info = PetscOptionsGetInt(PETSC_NULL,"-bmx",&this->bmx,&flg); CHKERRQ(info);
  info = PetscOptionsGetInt(PETSC_NULL,"-bmy",&this->bmy,&flg); CHKERRQ(info);
  //  info = PetscOptionsGetDouble(PETSC_NULL,"-bheight",&this->bheight,&flg); CHKERRQ(info);
  this->bmy = PetscMax(0,this->bmy);this->bmy = PetscMin(my,this->bmy);
  this->bmx = PetscMax(0,this->bmx);this->bmx = PetscMin(mx,this->bmx);
  bmy=this->bmy, bmx=this->bmx;

  info = DAGetCorners(this->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);

  info = XL->put(-20.0); CHKERRQ(info);
  info = XU->put(20.0); CHKERRQ(info);

  info = XL->getCoefPtrReadWriteLock(xl); CHKERRQ(info);

  info = PetscOptionsHasName(PETSC_NULL,"-cylinder",&cylinder); CHKERRQ(info);
  /* Compute the optional lower box */
  if (cylinder){
    for (i=xs; i< xs+xm; i++){    
      for (j=ys; j<ys+ym; j++){
        row=(j-ys)*xm + (i-xs);
        t1=(2.0*i-mx)*bmy;
        t2=(2.0*j-my)*bmx;
        t3=bmx*bmx*bmy*bmy;
        if ( t1*t1 + t2*t2 <= t3 ){
          xl[row] = this->bheight;
        }
      }
    }
  } else {
    for (i=xs; i< xs+xm; i++){    
      for (j=ys; j<ys+ym; j++){
        row=(j-ys)*xm + (i-xs);
        if (i>=(mx-bmx)/2 && i<mx-(mx-bmx)/2 && 
            j>=(my-bmy)/2 && j<my-(my-bmy)/2 ){
          xl[row] = this->bheight;
        }
      }
    }
  }

  info = XL->releaseCoefPtrLock(xl); CHKERRQ(info);
  return 0;
}

int MinsurfExample::getVariableVector(ESIVCTR **xx){
  *xx=this->xvec;
  return 0;
}

int MinsurfExample::getHessianMatrix(esi::Operator<double,int> **HH){
  *HH=this->H;
  return 0;
}
int MinsurfExample::monitor(){
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "MinsurfExample::evaluateObjectiveAndGradient"
/*  
    ESIFG - Evaluates the function, f(X), and gradient, G(X). 

    Input Parameters:
.   X   - input vector
    
    Output Parameters:
+   G - vector containing the newly evaluated gradient
-   f - function value

*/
int MinsurfExample::evaluateObjectiveAndGradient(ESIVCTR *XESI, double *fcn, ESIVCTR *GESI){

  int nn;

  int    info,i,j,row,row2;
  int    xs,xm,gxs,gxm,ys,ym,gys,gym;
  double zero=0.0,ft=0;
  double hx=1.0/(mx+1),hy=1.0/(my+1), hydhx=hy/hx, hxdhy=hx/hy, area=0.5*hx*hy;
  double rhx=mx+1, rhy=my+1;
  double f1,f2,f3,f4,f5,f6,d1,d2,d3,d4,d5,d6,d7,d8,xc,xl,xr,xt,xb,xlt,xrb;
  double df1dxc,df2dxc,df3dxc,df4dxc,df5dxc,df6dxc;
  double *g, *x;
  double *gesi,*xesi,*xglobal,*gglobal;

  /* Get local mesh boundaries */
  info = DAGetCorners(this->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(this->da,&gxs,&gys,PETSC_NULL,&gxm,&gym,PETSC_NULL); CHKERRQ(info);

  info = XESI->getCoefPtrReadWriteLock(xesi); CHKERRQ(info);
  info = VecGetLocalSize(globalX,&nn); CHKERRQ(info);
  info = VecGetArray(globalX,&xglobal); CHKERRQ(info);
  for (i=0;i<nn;i++){  xglobal[i]=xesi[i]; }
  info = VecRestoreArray(globalX,&xglobal); CHKERRQ(info);
  info = XESI->releaseCoefPtrLock(xesi); CHKERRQ(info);

  /* Scatter ghost points to local vector */
  info = DAGlobalToLocalBegin(this->da,globalX,INSERT_VALUES,localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(this->da,globalX,INSERT_VALUES,localX); CHKERRQ(info);

  /* Initialize vector to zero */
  info = VecSet(globalG, zero); CHKERRQ(info);

  /* Get pointers to vector data */
  info = VecGetArray(localX,&x); CHKERRQ(info);
  info = VecGetArray(globalG,&g); CHKERRQ(info);

  /* Compute function and gradient over the locally owned part of the mesh */
  for (j=ys; j<ys+ym; j++){
    for (i=xs; i< xs+xm; i++){
      row2=(j-ys)*xm + (i-xs);
      row=(j-gys)*gxm + (i-gxs);

      xc = x[row];
      xlt=xrb=xl=xr=xb=xt=xc;
      
      if (i==0){ /* left side */
        xl= this->left[j-ys+1];
        xlt = this->left[j-ys+2];
      } else {
        xl = x[row-1];
      }

      if (j==0){ /* bottom side */
        xb=this->bottom[i-xs+1];
        xrb =this->bottom[i-xs+2];
      } else {
        xb = x[row-gxm];
      }
      
      if (i+1 == gxs+gxm){ /* right side */
        xr=this->right[j-ys+1];
        xrb = this->right[j-ys];
      } else {
        xr = x[row+1];
      }

      if (j+1==gys+gym){ /* top side */
        xt=this->top[i-xs+1];
        xlt = this->top[i-xs];
      } else {
        xt = x[row+gxm];
      }

      if (i>gxs && j+1<gys+gym){
        xlt = x[row-1+gxm];
      }
      if (j>gys && i+1<gxs+gxm){
        xrb = x[row+1-gxm];
      }

      d1 = (xc-xl);
      d2 = (xc-xr);
      d3 = (xc-xt);
      d4 = (xc-xb);
      d5 = (xr-xrb);
      d6 = (xrb-xb);
      d7 = (xlt-xl);
      d8 = (xt-xlt);
      
      df1dxc = d1*hydhx;
      df2dxc = ( d1*hydhx + d4*hxdhy );
      df3dxc = d3*hxdhy;
      df4dxc = ( d2*hydhx + d3*hxdhy );
      df5dxc = d2*hydhx;
      df6dxc = d4*hxdhy;

      d1 *= rhx;
      d2 *= rhx;
      d3 *= rhy;
      d4 *= rhy;
      d5 *= rhy;
      d6 *= rhx;
      d7 *= rhy;
      d8 *= rhx;

      f1 = sqrt( 1.0 + d1*d1 + d7*d7);
      f2 = sqrt( 1.0 + d1*d1 + d4*d4);
      f3 = sqrt( 1.0 + d3*d3 + d8*d8);
      f4 = sqrt( 1.0 + d3*d3 + d2*d2);
      f5 = sqrt( 1.0 + d2*d2 + d5*d5);
      f6 = sqrt( 1.0 + d4*d4 + d6*d6);
      
      f2 = sqrt( 1.0 + d1*d1 + d4*d4);
      f4 = sqrt( 1.0 + d3*d3 + d2*d2);

      ft = ft + (f2 + f4);

      df1dxc /= f1;
      df2dxc /= f2;
      df3dxc /= f3;
      df4dxc /= f4;
      df5dxc /= f5;
      df6dxc /= f6;

      g[row2] = (df1dxc+df2dxc+df3dxc+df4dxc+df5dxc+df6dxc ) * 0.5;

    }
  }

  /* Compute triangular areas along the border of the domain. */
  if (xs==0){ /* left side */
    for (j=ys; j<ys+ym; j++){
      d3=(this->left[j-ys+1] - this->left[j-ys+2])*rhy;
      d2=(this->left[j-ys+1] - x[(j-gys)*gxm])*rhx;
      ft = ft+sqrt( 1.0 + d3*d3 + d2*d2);
    }
  }

  if (ys==0){ /* bottom side */
    for (i=xs; i<xs+xm; i++){
      d2=(this->bottom[i+1-xs]-this->bottom[i-xs+2])*rhx;
      d3=(this->bottom[i-xs+1]-x[i-gxs])*rhy;
      ft = ft+sqrt( 1.0 + d3*d3 + d2*d2);
    }
  }

  if (xs+xm==mx){ /* right side */
    for (j=ys; j< ys+ym; j++){
      d1=(x[(j+1-gys)*gxm-1]-this->right[j-ys+1])*rhx;
      d4=(this->right[j-ys]-this->right[j-ys+1])*rhy;
      ft = ft+sqrt( 1.0 + d1*d1 + d4*d4);
    }
  }
  if (ys+ym==my){ /* top side */
    for (i=xs; i<xs+xm; i++){
      d1=(x[(gym-1)*gxm + i-gxs] - this->top[i-xs+1])*rhy;
      d4=(this->top[i-xs+1] - this->top[i-xs])*rhx;
      ft = ft+sqrt( 1.0 + d1*d1 + d4*d4);
    }
  }

 if (ys==0 && xs==0){
    d1=(this->left[0]-this->left[1])/hy;
    d2=(this->bottom[0]-this->bottom[1])*rhx;
    ft +=sqrt( 1.0 + d1*d1 + d2*d2);
  }
  if (ys+ym == my && xs+xm == mx){
    d1=(this->right[ym+1] - this->right[ym])*rhy;
    d2=(this->top[xm+1] - this->top[xm])*rhx;
    ft +=sqrt( 1.0 + d1*d1 + d2*d2);
  }

  ft=ft*area;
  info = MPI_Allreduce(&ft,fcn,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);CHKERRQ(info);

  /* Restore vectors */
  info = VecRestoreArray(localX,&x); CHKERRQ(info);
  info = VecRestoreArray(globalG,&g); CHKERRQ(info);

  info = VecGetLocalSize(globalG,&nn); CHKERRQ(info);
  info = GESI->getCoefPtrReadWriteLock(gesi); CHKERRQ(info);
  info = VecGetArray(globalG,&gglobal); CHKERRQ(info);
  for (i=0;i<nn;i++){ gesi[i] = gglobal[i]; /* printf("%d  %4.4f, ", i,gesi[i]); */ }
  info = VecRestoreArray(globalG,&gglobal); CHKERRQ(info);
  info = GESI->releaseCoefPtrLock(gesi); CHKERRQ(info);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "MinsurfExample::evaluateHessian"
int MinsurfExample::evaluateHessian(ESIVCTR *XESI, ESIMAT *HH){

  esi::MatrixRowWriteAccess<double,int>* HHH;
  esi::IndexSpace<int>* rowmap = NULL;
  int    i,j,k,row,info;
  int    xs,xm,gxs,gxm,ys,ym,gys,gym,col[7];
  double hx=1.0/(mx+1), hy=1.0/(my+1), hydhx=hy/hx, hxdhy=hx/hy;
  double f1,f2,f3,f4,f5,f6,d1,d2,d3,d4,d5,d6,d7,d8,xc,xl,xr,xt,xb,xlt,xrb;
  double hl,hr,ht,hb,hc,htl,hbr;
  double *x, v[7];
  double *xesi;
  int nn;
  double *xglobal, *globalrow;
  int firstLocalEqn;
  
  /* Get local mesh boundaries */
  info = DAGetCorners(this->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(this->da,&gxs,&gys,PETSC_NULL,&gxm,&gym,PETSC_NULL); CHKERRQ(info);

  //  esi::MatrixRowWriteAccess<double,int>* HHH;
  HHH = dynamic_cast<esi::MatrixRowWriteAccess<double,int>*>(HH);

  info = XESI->getIndexSpace(rowmap); CHKERRQ(info);
  info = rowmap->getLocalPartitionOffset(firstLocalEqn); CHKERRQ(info);

  info = XESI->getCoefPtrReadWriteLock(xesi); CHKERRQ(info);
  info = VecGetLocalSize(globalX,&nn); CHKERRQ(info);
  info = VecGetArray(globalX,&xglobal); CHKERRQ(info);
  for (i=0;i<nn;i++){  xglobal[i]=xesi[i]; }
  info = VecRestoreArray(globalX,&xglobal); CHKERRQ(info);
  info = XESI->releaseCoefPtrLock(xesi); CHKERRQ(info);

  /* Scatter ghost points to local vector */
  info = DAGlobalToLocalBegin(this->da,globalX,INSERT_VALUES,localX); CHKERRQ(info);
  info = DAGlobalToLocalEnd(this->da,globalX,INSERT_VALUES,localX); CHKERRQ(info);

  /* Get pointers to vector data */
  info = VecGetArray(orderX,&globalrow); CHKERRQ(info);
  info = VecGetArray(localX,&x); CHKERRQ(info);

  /* Compute Hessian over the locally owned part of the mesh */

  for (j=0; j<ym; j++){
    for (i=0; i<xm; i++){

      row=(j+ys-gys)*gxm + (i+xs-gxs);
      xc = x[row]; 
      xlt=xrb=xl=xr=xb=xt=xc;

      /* Left side */
      if (i+xs==gxs){
        xl= this->left[j+1];
        xlt = this->left[j+2];
      } else {
        xl = x[row-1];
      }
      
      if (j+ys==gys){
        xb=this->bottom[i+1];
        xrb = this->bottom[i+2];
      } else {
        xb = x[row-gxm];
      }
      
      if (i+xs+1 == gxs+gxm){
        xr=this->right[j+1];
        xrb = this->right[j];
      } else {
        xr = x[row+1];
      }

      if (j+ys+1==gys+gym){
        xt=this->top[i+1];
        xlt = this->top[i];
      }else {
        xt = x[row+gxm];
      }

      if (i+xs>gxs && j+ys+1<gys+gym){
        xlt = x[row-1+gxm];
      }
      if (j+ys>gys && i+xs+1<gxs+gxm){
        xrb = x[row+1-gxm];
      }


      d1 = (xc-xl)/hx;
      d2 = (xc-xr)/hx;
      d3 = (xc-xt)/hy;
      d4 = (xc-xb)/hy;
      d5 = (xrb-xr)/hy;
      d6 = (xrb-xb)/hx;
      d7 = (xlt-xl)/hy;
      d8 = (xlt-xt)/hx;
      
      f1 = sqrt( 1.0 + d1*d1 + d7*d7);
      f2 = sqrt( 1.0 + d1*d1 + d4*d4);
      f3 = sqrt( 1.0 + d3*d3 + d8*d8);
      f4 = sqrt( 1.0 + d3*d3 + d2*d2);
      f5 = sqrt( 1.0 + d2*d2 + d5*d5);
      f6 = sqrt( 1.0 + d4*d4 + d6*d6);


      hl = (-hydhx*(1.0+d7*d7)+d1*d7)/(f1*f1*f1)+
	(-hydhx*(1.0+d4*d4)+d1*d4)/(f2*f2*f2);
      hr = (-hydhx*(1.0+d5*d5)+d2*d5)/(f5*f5*f5)+
	(-hydhx*(1.0+d3*d3)+d2*d3)/(f4*f4*f4);
      ht = (-hxdhy*(1.0+d8*d8)+d3*d8)/(f3*f3*f3)+
	(-hxdhy*(1.0+d2*d2)+d2*d3)/(f4*f4*f4);
      hb = (-hxdhy*(1.0+d6*d6)+d4*d6)/(f6*f6*f6)+
	(-hxdhy*(1.0+d1*d1)+d1*d4)/(f2*f2*f2);

      hbr = -d2*d5/(f5*f5*f5) - d4*d6/(f6*f6*f6);
      htl = -d1*d7/(f1*f1*f1) - d3*d8/(f3*f3*f3);

      hc = hydhx*(1.0+d7*d7)/(f1*f1*f1) + hxdhy*(1.0+d8*d8)/(f3*f3*f3) +
	hydhx*(1.0+d5*d5)/(f5*f5*f5) + hxdhy*(1.0+d6*d6)/(f6*f6*f6) +
	(hxdhy*(1.0+d1*d1)+hydhx*(1.0+d4*d4)-2*d1*d4)/(f2*f2*f2) +
	(hxdhy*(1.0+d2*d2)+hydhx*(1.0+d3*d3)-2*d2*d3)/(f4*f4*f4);

      hl/=2.0; hr/=2.0; ht/=2.0; hb/=2.0; hbr/=2.0; htl/=2.0;  hc/=2.0; 

      k=0;
      col[0]=-1;col[1]=-1;col[2]=-1;col[3]=-1;col[4]=-1;col[5]=-1;col[6]=-1;
      v[0]=0.0; v[1]=0.0; v[2]=0.0; v[3]=0.0; v[4]=0.0; v[5]=0.0; v[6]=0.0;
      if (j+ys>0){ 
	v[k]=hb; col[k]=(int)globalrow[row-xm]; k++;
      }
      
      if (j+ys>0 && i+xs < mx -1){
	v[k]=hbr; col[k]=(int)globalrow[row-xm+1]; k++;
      }
      
      if (i+xs>0){
	v[k]= hl; col[k]=(int)globalrow[row-1]; k++;
      }
      
      v[k]= hc; col[k]=(int)globalrow[row]; k++;
      
      if (i+xs < mx-1 ){
	v[k]= hr; col[k]=(int)globalrow[row+1]; k++;
      }
      
      if (i+xs>0 && j+ys < my-1 ){
	v[k]= htl; col[k] = (int)globalrow[row+xm-1]; k++;
      }
      
      if (j+ys < my-1 ){
	v[k]= ht; col[k] =(int)globalrow[row+xm]; k++;
      }

      /* 
	 Set matrix values using global numbering
      */
      info = HHH->copyIntoRow(globalrow[row], v, col, k); CHKERRQ(info);
      
    }
  }
   
  /* Restore vectors */
  info = VecRestoreArray(localX,&x); CHKERRQ(info);
  info = VecRestoreArray(orderX,&globalrow); CHKERRQ(info);

  CHK_ERR( HH->setup() );


  //  esi::MatrixRowReadAccess<double,int>* HR;
  //  HR = dynamic_cast<esi::MatrixRowReadAccess<double,int>*>(HH);
  //  info = cout_ESI_RowMatrix(HR);
 
  return 0;
}



#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  TAO_SOLVER tao;                   /* TAO_SOLVER solver context */
  int        info;                  /* error code */
  TaoTerminateReason reason;
  MinsurfExample *myesiapp;
  TaoESIApplication *mytaoapp;
  esi::Solver<double,int>* lsolver=NULL;

  /* Initialize TAO */
  PetscInitialize(&argc,&argv,(char *)0,help);
  TaoInitialize(&argc,&argv,(char *)0,help);

  myesiapp = new MinsurfExample();
  mytaoapp = new TaoESIApplication(myesiapp);

  /* Create TAO solver */
  info = TaoCreate(MPI_COMM_WORLD,"tao_lmvm",&tao); CHKERRQ(info);

  TaoLinearSolverESI *tsolver=NULL;
  info = create_ESI_Solver(myesiapp->H, lsolver);
  info = TaoWrapESIKSP(lsolver, &tsolver); CHKERRQ(info);
  info = TaoSetLinearSolver(tao,tsolver);CHKERRQ(info);

  /* Set routine for function and gradient evaluation */
  info = TaoSetApplication(tao,mytaoapp);  CHKERRQ(info);


  info = TaoSetOptions(tao);CHKERRQ(info);
  // info = TaoSetInitialVector(tao,xx); CHKERRQ(info);
  info = TaoSolve(tao); CHKERRQ(info);

  info = TaoGetIterationData(tao,0,0,0,0,0,&reason); CHKERRQ(info);
  if (reason <= 0)
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");

  /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
 
  delete mytaoapp;
  delete myesiapp;
  delete lsolver;

  /* Finalize TAO */
  TaoFinalize();
  PetscFinalize();
  return 0;
}

//----------------------------------------------
int create_ESI_IndexSpace(int numGlobal, int numLocal,
                   Epetra_Comm& comm, esi::IndexSpace<int>*& esi_indexspace)
{
   //we're using indexBase = 0...
   esi_indexspace = new epetra_esi::IndexSpace<int>(numGlobal, numLocal, 0, comm);

   //return error-code -1 if the allocation failed
   if (esi_indexspace == NULL) return(-1);

   return(0);
}

//----------------------------------------------
int create_ESI_Vector(esi::IndexSpace<int>* indexspace,
                      esi::Vector<double,int>*& vec)
{

  //In order to construct a epetra_esi::Vector, we need an actual
  //epetra_esi::IndexSpace...
  //
  epetra_esi::IndexSpace<int>* petraindexspace = NULL;
  CHK_ERR( indexspace->getInterface("epetra_esi::IndexSpace", (void*&)petraindexspace) );

  vec = new epetra_esi::Vector<double,int>(*petraindexspace);
  if (vec == NULL) return(-1);
  return(0);
}


/* -------------------------------------------------------------------- */
int  MinsurfExample::create_MESI_Operator(int numGlobal, int numLocal,
					  esi::IndexSpace<int>* indexspace)
{
  //The goal of this function is to create an esi::Operator. But since the
  //run-time type of the esi::Operator will be Petra_ESI_CRS_Matrix, we first
  //have to create and fill a Petra_CRS_Graph (the container that defines the
  //structure of the matrix.
  //
  int i,j,k,row,info,err;
  int    xs,xm,gxs,gxm,ys,ym,gys,gym,col[7];
  double *globalrow;
  double v[7];
  Epetra_CrsGraph* graph;

  //We'll start by finding out what the first local equation-number is.
  //
  /* Get local mesh boundaries */
  info = DAGetCorners(this->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(this->da,&gxs,&gys,PETSC_NULL,&gxm,&gym,PETSC_NULL); CHKERRQ(info);

  int localOffset = 0;
  CHK_ERR( indexspace->getLocalPartitionOffset(localOffset) );

  //This will be a very simple structure.
  //We want our matrix to have a diagonal, and the first off-diagonals both
  //above and below the main diagonal.
  //
  int* rowLengths = new int[numLocal];
  if (rowLengths == NULL) return(-1);

  for(i=0; i<numLocal; i++) rowLengths[i] = 7;
  if (localOffset == 0) rowLengths[0] = 7;
  if (localOffset+numLocal == numGlobal) rowLengths[numLocal-1] = 7;

  //We need a native Petra_Map to construct the Petra_CRS_Graph with.
  Epetra_Map* petramap = NULL;
  CHK_ERR( indexspace->getInterface("Epetra_Map", (void*&)petramap) );

  graph = new Epetra_CrsGraph(Copy, *petramap, rowLengths);

  delete[] rowLengths;
  info = VecGetArray(orderX,&globalrow); CHKERRQ(info);

  //Now we've allocated the 'shell' of the graph, so let's proceed with
  //setting the column-indices...


  for (j=0; j<ym; j++){

    for (i=0; i< xm; i++){

      row=(j+ys-gys)*gxm + (i+xs-gxs);
      k=0;

      if (j+ys>0){ 
	col[k]=(int)globalrow[row - gxm]; k++;
      }
      
      if (j+ys>0 &&  i+xs < mx -1){
	col[k]=(int)globalrow[row - gxm+1]; k++;
      }
      
      if (i+xs>0){
	col[k]=(int)globalrow[row - 1]; k++;
      }
      
      v[k]= col[k]=(int)globalrow[row]; k++;
      
      if (i+xs < mx-1 ){
	col[k]=(int)globalrow[row+1]; k++;
      }
      
      if (j+ys < my-1 ){
	col[k] = (int)globalrow[row+gxm]; k++;
      }

      if (i+xs>0 && j+ys < my-1 ){
	col[k] = (int)globalrow[row+gxm-1]; k++;
      }
      
      err = graph->InsertGlobalIndices(globalrow[row], k, col);
      if (err != 0) return(err);

    }
  }
  info = VecRestoreArray(orderX,&globalrow); CHKERRQ(info);

  err = graph->TransformToLocal();
  if (err != 0) return(-1);

  //Now we've got a fully initialized graph holding the structure of our
  //matrix, so we're ready to go ahead and construct the matrix.
  //
  epetra_esi::CrsMatrix<double,int>* petramat =   new epetra_esi::CrsMatrix<double,int>(Copy, *graph);
  if (petramat == NULL) return(-1);
  this->petraH=petramat;

  CHK_ERR( petramat->getInterface("esi::Operator", (void*&)this->H) );
  CHK_ERR( this->H->setup() );

  return(0);
}




int create_ESI_Solver(esi::Operator<double,int>* A,
		      esi::Solver<double,int>*& solver)
{
  //The run-time type of the esi::Solver will be Aztec_ESI_Solver. We
  //
  epetra_esi::CrsMatrix<double,int>* petraA = NULL;
  CHK_ERR( A->getInterface("epetra_esi::CrsMatrix", (void*&)petraA) );

  solver = new aztecoo_esi::Solver<double,int>(petraA);
  if (solver == NULL) return(-1);

  /*
  return 0;
  solver = new aztecoo_esi::Solver<double,int>(A);
  if (solver == NULL) return(-1);
  */
   int numParams = 8;
   char** paramStrings = new char*[numParams];
   for(int i=0; i<numParams; i++) paramStrings[i] = new char[128];
   sprintf(paramStrings[0], "AZ_solver AZ_gmres");
   sprintf(paramStrings[1], "AZ_tol 1.e-7");
   sprintf(paramStrings[2], "AZ_max_iter 200");
   sprintf(paramStrings[3], "AZ_output AZ_all");
//   sprintf(paramStrings[4], "AZ_precond AZ_Neumann");
   sprintf(paramStrings[4], "AZ_precond AZ_dom_decomp");
   sprintf(paramStrings[5], "AZ_subdomain_solver AZ_ilut");
   sprintf(paramStrings[6], "AZ_ilut_fill 0.8");
   sprintf(paramStrings[7], "AZ_drop 0.01");

   CHK_ERR( solver->parameters(numParams, paramStrings) );
   for(int j=0; j<numParams; j++) delete [] paramStrings[j];
   delete [] paramStrings;

  return(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MSA_BoundaryConditions"
/* 
   MSA_BoundaryConditions -  Calculates the boundary conditions for
   the region.

   Input Parameter:
.  user - user-defined application context

   Output Parameter:
.  user - user-defined application context
*/
int MinsurfExample::MSA_BoundaryConditions()
{
  int        i,j,k,limit=0,info,maxits=5;
  int        xs,ys,xm,ym,gxs,gys,gxm,gym;
  int        mx=this->mx,my=this->my;
  int        bsize=0, lsize=0, tsize=0, rsize=0;
  double     one=1.0, two=2.0, three=3.0, tol=1e-10;
  double     fnorm,det,hx,hy,xt=0,yt=0;
  double     u1,u2,nf1,nf2,njac11,njac12,njac21,njac22;
  double     b=-0.5, t=0.5, l=-0.5, r=0.5;
  double     *boundary;
  PetscTruth   flg;

  /* Get local mesh boundaries */
  info = DAGetCorners(this->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL); CHKERRQ(info);
  info = DAGetGhostCorners(this->da,&gxs,&gys,PETSC_NULL,&gxm,&gym,PETSC_NULL); CHKERRQ(info);

  bsize=xm+2;
  lsize=ym+2;
  rsize=ym+2;
  tsize=xm+2;

  info = PetscMalloc(bsize*sizeof(double),&this->bottom); CHKERRQ(info);
  info = PetscMalloc(tsize*sizeof(double),&this->top); CHKERRQ(info);
  info = PetscMalloc(lsize*sizeof(double),&this->left); CHKERRQ(info);
  info = PetscMalloc(rsize*sizeof(double),&this->right); CHKERRQ(info);

  hx= (r-l)/(mx+1); hy=(t-b)/(my+1);

  for (j=0; j<4; j++){
    if (j==0){
      yt=b;
      xt=l+hx*xs;
      limit=bsize;
      boundary=this->bottom;
    } else if (j==1){
      yt=t;
      xt=l+hx*xs;
      limit=tsize;
      boundary=this->top;
    } else if (j==2){
      yt=b+hy*ys;
      xt=l;
      limit=lsize;
      boundary=this->left;
    } else if (j==3){
      yt=b+hy*ys;
      xt=r;
      limit=rsize;
      boundary=this->right;
    }

    for (i=0; i<limit; i++){
      u1=xt;
      u2=-yt;
      for (k=0; k<maxits; k++){
	nf1=u1 + u1*u2*u2 - u1*u1*u1/three-xt;
	nf2=-u2 - u1*u1*u2 + u2*u2*u2/three-yt;
	fnorm=sqrt(nf1*nf1+nf2*nf2);
	if (fnorm <= tol) break;
	njac11=one+u2*u2-u1*u1;
	njac12=two*u1*u2;
	njac21=-two*u1*u2;
	njac22=-one - u1*u1 + u2*u2;
	det = njac11*njac22-njac21*njac12;
	u1 = u1-(njac22*nf1-njac12*nf2)/det;
	u2 = u2-(njac11*nf2-njac21*nf1)/det;
      }

      boundary[i]=u1*u1-u2*u2;
      if (j==0 || j==1) {
	xt=xt+hx;
      } else if (j==2 || j==3){
	yt=yt+hy;
      }
      
    }

  }

  /* Scale the boundary if desired */
  if (1==1){ 
    double scl = 1.0;

    info = PetscOptionsGetReal(PETSC_NULL,"-bottom",&scl,&flg); 
    CHKERRQ(info);
    if (flg){
      for (i=0;i<bsize;i++) this->bottom[i]*=scl;
    }

    info = PetscOptionsGetReal(PETSC_NULL,"-top",&scl,&flg); 
    CHKERRQ(info);
    if (flg){
      for (i=0;i<tsize;i++) this->top[i]*=scl;
    }

    info = PetscOptionsGetReal(PETSC_NULL,"-right",&scl,&flg); 
    CHKERRQ(info);
    if (flg){
      for (i=0;i<rsize;i++) this->right[i]*=scl;
    }

    info = PetscOptionsGetReal(PETSC_NULL,"-left",&scl,&flg); 
    CHKERRQ(info);
    if (flg){
      for (i=0;i<lsize;i++) this->left[i]*=scl;
    }
  }
  
  return 0;
}

//----------------------------------------------
int cout_ESI_RowMatrix(esi::MatrixRowReadAccess<double,int>* esimat)
{
   esi::IndexSpace<int>* rowIS = NULL, *colIS = NULL;
   CHK_ERR(esimat->getIndexSpaces(rowIS, colIS) );

   //we need to know how many local equations there are, and what the first
   //local equation is. (Since we happen to know that the index-spaces were
   //built with an indexBase of 0,
   // then firstLocalEqn == 'getLocalPartitionOffset'.)
   int localSize, firstLocalEqn, localRank;
   CHK_ERR( rowIS->getLocalSize(localSize) );
   CHK_ERR( rowIS->getLocalPartitionOffset(firstLocalEqn) );
   CHK_ERR( rowIS->getLocalPartitionRank(localRank) );

   int numRows, numCols;
   CHK_ERR( esimat->getGlobalSizes(numRows, numCols) );
   cout << localRank << ": global rows " << numRows << ", global cols "
     << numCols << endl;

   Epetra_Array<int> colIndices;
   Epetra_Array<double> coefs;

   for(int i=0; i<localSize; i++) {
      //first, make sure our colIndices and coefs arrays are the right length.
      //(This operation is a no-op if the array is already the right size or
      //bigger.)
      int rowLen;
      int row = firstLocalEqn+i;
      CHK_ERR( esimat->getRowNonzeros(row, rowLen) );

      CHK_ERR( colIndices.resize(rowLen) );
      CHK_ERR( coefs.resize(rowLen) );

      CHK_ERR( esimat->copyOutRow(row, coefs.dataPtr(), colIndices.dataPtr(),
                                  coefs.length(), rowLen) );
      cout << localRank << ": row " << row << ": ";
      for(int j=0; j<rowLen; j++) {
         cout << "("<<colIndices[j]<<","<<coefs[j]<<") ";
      }
      cout << endl;
   }

   return(0);
}

//----------------------------------------------
int cout_ESI_Vector(esi::Vector<double,int>* esivec)
{
   esi::IndexSpace<int>* indexspace = NULL;
   CHK_ERR(esivec->getIndexSpace(indexspace) );

   //we need to know how many local equations there are, and what the first
   //local equation is. (Since we happen to know that the index-space was built
   //with an indexBase of 0, then firstLocalEqn == 'getLocalPartitionOffset'.)
   int localSize, firstLocalEqn, localRank;
   CHK_ERR( indexspace->getLocalSize(localSize) );
   CHK_ERR( indexspace->getLocalPartitionOffset(firstLocalEqn) );
   CHK_ERR( indexspace->getLocalPartitionRank(localRank) );

   double* coefs;
   CHK_ERR( esivec->getCoefPtrReadLock(coefs) );

   for(int i=0; i<localSize; i++) {
      cout << localRank << ": "<<firstLocalEqn+i<<", "<<coefs[i]<<endl;
   }

   CHK_ERR( esivec->releaseCoefPtrLock(coefs) );

   return(0);
}

