/* Program usage: coating [-help] [all TAO options] */

/* 
   Include "tao.h" so that we can use TAO solvers.  Note that this
   file automatically includes libraries such as:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - sysem routines        petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

*/

#include "tao.h"


static char help[] = 
"This example demonstrates use of the TAO package to solve a nonlinear \n\
least squares problem.  This example is based on a problem from the\n\
MINPACK-2 test suit.  The Chebyshev problem arises from the determination\n\
of the nodes of a quadature formula with equal weights. \n\n";


/* T
   Concepts: TAO - Solving a system of nonlinear equations, nonlinear least squares;
   Routines: TaoInitialize(); TaoFinalize(); 
   Routines: TaoCreate();  TaoDestroy();
   Routines: TaoPetscApplicationCreate(); TaoApplicationDestroy();
   Routines: TaoSetPetscFunction(); TaoSetPetscJacobian(); 
   Routines: TaoSetPetscConstraintsFunction();
   Routines: TaoSetPetscInitialVector(); TaoSetApplication();
   Routines: TaoSolve(); TaoSetFromOptions(); 
   Processors: 1
T*/ 



/* User-defined application context */
typedef struct {
  /* Application parameters */
  int nvar;     /* number of variables */
  int ncnst;    /* number of constraints */
  int mdiv2, mdiv4;
  double scale1, scale2;

  double *y, *V0, *V1; /* application data */
    
} AppCtx;

/* User provided Routines */
int InitializeVectors(Vec v, Vec Y);
int FormStartingPoint(AppCtx *user, Vec X);
int EvaluateConstraints(TAO_SOLVER tao, Vec Y, Vec fvec, void *user);
int EvaluateJacobian(TAO_SOLVER tao, Vec X, Mat *J, void *ptr);
int AppCtxInitialize(AppCtx*);
int AppCtxDestroy(AppCtx*);

/*--------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  int        info;              /* used to check for functions returning nonzeros */
  Vec        x, f;              /* solution, function */
  Mat        J;                 /* Jacobian of constraints function */
  TAO_SOLVER tao;               /* TAO_SOLVER solver context */
  TAO_APPLICATION taoapp;       /* the PETSc application */
  int*    matarray;
  TaoMethod  method = "tao_nlsq"; /* minimization method */
  int        i;
  AppCtx     user;              /* application context */

  /* Initialize TAO and PETSc */
  PetscInitialize(&argc,&argv,(char *)0,help);
  TaoInitialize(&argc,&argv,(char *)0,help);

  /* Initialize application data -- user function */
  info = AppCtxInitialize(&user); CHKERRQ(info);

  /* Allocate vectors */
  info = VecCreateSeq(MPI_COMM_SELF,user.nvar,&x); CHKERRQ(info);
  info = VecCreateSeq(MPI_COMM_SELF,user.ncnst,&f); CHKERRQ(info);

  /* Create Jacobian matrix  */
  info = PetscMalloc(sizeof(int)*user.ncnst, (void*) &matarray);CHKERRQ(info);
  for(i=0; i< user.mdiv2;i++){
    matarray[i]=6; matarray[i+user.mdiv2]=1;
  }
    
  info = MatCreateSeqAIJ(MPI_COMM_WORLD,user.ncnst,user.nvar,TAO_NULL,matarray,&J);
  CHKERRQ(info);

  /* TAO code begins here */

  /* Create TAO solver */
  info = TaoCreate(MPI_COMM_SELF,method,&tao);CHKERRQ(info);
  info = TaoPetscApplicationCreate(PETSC_COMM_SELF,&taoapp); CHKERRQ(info);


  /* Set the function and Jacobian routines. */
  info = TaoSetPetscFunction(taoapp,x,TAO_NULL,TAO_NULL); 
  CHKERRQ(info);
  info = TaoSetPetscJacobian(taoapp, J, EvaluateJacobian, (void*)&user); 
  CHKERRQ(info);
  info = TaoSetPetscConstraintsFunction(taoapp, f, EvaluateConstraints,(void*) &user); 
  CHKERRQ(info);

   /*  Compute the standard starting point. */
  info = FormStartingPoint(&user, x); CHKERRQ(info);
  info = TaoSetPetscInitialVector(taoapp,x); CHKERRQ(info);

  /* Now that the PETSc application is set, attach to TAO context */
  info = TaoSetApplication(tao,taoapp); CHKERRQ(info); 

  /* Check for any TAO command line arguments */
  info = TaoSetFromOptions(tao); CHKERRQ(info);

  /* Perform the Solve */
  info = TaoSolve(tao); CHKERRQ(info);


  /* Free TAO data structures */
  info = TaoApplicationDestroy(taoapp); CHKERRQ(info);
  info = TaoDestroy(tao); CHKERRQ(info);

   /* Free PETSc data structures */
  info = VecDestroy(x); CHKERRQ(info);
  info = VecDestroy(f); CHKERRQ(info);
  info = MatDestroy(J); CHKERRQ(info);

  /* Free user data structures */
  info = AppCtxDestroy(&user); CHKERRQ(info);
  info = PetscFree(matarray);CHKERRQ(info);


  /* Finalize TAO */
  PetscFinalize();
  TaoFinalize();

  return 0;     
}


/*--------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "EvaluateConstraints"
int EvaluateConstraints(TAO_SOLVER tao, Vec X, Vec F, void *ptr){

  AppCtx *user = (AppCtx *)ptr;
  int i,info;
  double zero=0.0,*xvec, *fvec;

  info = VecSet(F, zero); CHKERRQ(info);
  info = VecGetArray(X,&xvec); CHKERRQ(info);
  info = VecGetArray(F,&fvec); CHKERRQ(info);

  for(i=0;i<user->mdiv4;i++){

    fvec[i] = xvec[0]+xvec[1]*(user->V0[i]+xvec[8+i])+xvec[2]*
      (user->V1[i]+xvec[8+i+user->mdiv4])+xvec[3]*(user->V0[i]+xvec[8+i])*
      (user->V1[i]+xvec[8+i+user->mdiv4])-user->y[i];
    
    fvec[i+user->mdiv4] = xvec[4]+xvec[5]*(user->V0[i]+xvec[8+i])+xvec[6]*
      (user->V1[i]+xvec[8+i+user->mdiv4])+xvec[7]*(user->V0[i]+xvec[8+i])* 
      (user->V1[i]+xvec[8+i+user->mdiv4]) - user->y[i+user->mdiv4];
    
    fvec[i+2*user->mdiv4] = user->scale1*xvec[8+i];
    
    fvec[i+3*user->mdiv4] = user->scale2*xvec[8+i+user->mdiv4];
  }


  info = VecRestoreArray(X,&xvec); CHKERRQ(info);
  info = VecRestoreArray(F,&fvec); CHKERRQ(info);

  PetscLogFlops(user->mdiv4 * 26);
  return 0;
}



/*--------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "EvaluateJacobian"
int EvaluateJacobian(TAO_SOLVER tao, Vec X, Mat *J, void *ptr){

  AppCtx     *user = (AppCtx *)ptr;
  int        i, j, info;
  int        *row, *column;
  double     zero = 0.0, one = 1.0;
  double     *xvec, **jac;

  info = MatZeroEntries(*J); CHKERRQ(info);
  info = PetscMalloc(sizeof(double*)*user->ncnst, (void*) &jac ); CHKERRQ(info);

  info = PetscMalloc(sizeof(double)*user->ncnst*user->nvar,(void*) &jac[0]); CHKERRQ(info);
  for (i=0;i<user->ncnst*user->nvar;i++) jac[0][i]=0.0;

  for(i=1;i<user->ncnst;i++){
    jac[i] = jac[i-1] + user->nvar;
  }

  info = PetscMalloc(sizeof(int)*user->ncnst,(void*) &row); CHKERRQ(info);
  info = PetscMalloc(sizeof(int)*user->nvar,(void*) &column); CHKERRQ(info);

  for(j=0; j < user->ncnst; j++){
    row[j] = j;
  }
  for(j=0;j< user->nvar; j++){
    column[j] = j;
  }

 for(j=0;j<(8+user->mdiv2);j++){
    for(i=0;i<user->nvar;i++){
      jac[i][j] = zero;
    }		
  }

  info = VecGetArray(X,&xvec); CHKERRQ(info);

  for(i=0;i<user->mdiv4;i++){

    jac[i][0]=one;
    jac[i][1]=user->V0[i] + xvec[8+i];
    jac[i][2]=user->V1[i] + xvec[8+i+user->mdiv4];
    jac[i][3]=(user->V0[i]+xvec[8+i])*(user->V1[i]+xvec[8+i+user->mdiv4]);
    jac[i][8+i]=xvec[1] + xvec[3]*(user->V1[i]+xvec[8+i+user->mdiv4]);
    jac[i][8+i+user->mdiv4]=xvec[2] + xvec[3]*(user->V0[i]+xvec[8+i]);
    jac[i+user->mdiv4][4]=one;
    jac[i+user->mdiv4][5]=user->V0[i] + xvec[8+i];
    jac[i+user->mdiv4][6]=user->V1[i] + xvec[8+i+user->mdiv4];
    jac[i+user->mdiv4][7]=(user->V0[i]+xvec[8+i])*
      (user->V1[i]+ xvec[8+i+user->mdiv4]);
    jac[i+user->mdiv4][8+i]=xvec[5] + xvec[7]*
      (user->V1[i]+xvec[8+i+user->mdiv4]);
    jac[i+user->mdiv4][8+i+user->mdiv4]=xvec[6] + xvec[7]*
      (user->V0[i]+xvec[8+i]);
    jac[i+2*user->mdiv4][8+i]=user->scale1;
    jac[i+3*user->mdiv4][8+i+user->mdiv4]=user->scale2;
  }

  info = MatSetValues(*J,user->ncnst, row, user->nvar, column, *jac,INSERT_VALUES);  CHKERRQ(info); 

  /*  for(i=0;i<user->ncnst;i++){
      for(j=0;j<15;j++){ 
      PetscPrintf(MPI_COMM_SELF,"jac[%d][0] = %f", i, jac[i][0]);
      }
      PetscPrintf(MPI_COMM_SELF,"\n"); 
  }*/


  info = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
 
  /*info = MatView(*J,VIEWER_STDOUT_SELF); CHKERRQ(info);   */

  info = VecRestoreArray(X,&xvec); CHKERRQ(info);

  PetscLogFlops(user->mdiv4 * 22);

  PetscFree(jac[0]);
  PetscFree(jac);
  PetscFree(row);
  PetscFree(column);

  return 0;

}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormStartingPoint"
int FormStartingPoint(AppCtx *ptr,Vec X)
{
  AppCtx *user = (AppCtx *) ptr;
  int    i, info;
  double zero = 0.0;
  double *x;

  info = VecGetArray(X,&x); CHKERRQ(info);

  x[0] = -8.0e0;
  x[1] = 13.0e0;
  x[2] = 1.2e0;
  x[3] = 0.2e0;
  x[4] = 0.1e0;
  x[5] = 6.0e0;
  x[6] = 5.5e0;
  x[7] = -5.2e0;
  for (i=0; i< user->mdiv2; i++)   x[8+i] = zero;

  info = VecRestoreArray(X,&x); CHKERRQ(info);

  return 0;

}




/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "AppCtxInitialize"
/*
    AppCtxInitialize - Sets up user-defined work space; initializes problem
    parameters.

    Input Parameter:
    user - user-defined application context

    Output Parameter:
    user - user-defined application context
*/
int AppCtxInitialize(AppCtx *user)
{
  int      info;

  user->ncnst = 252; user->nvar = 134;
  user->mdiv2=126; user->mdiv4=63;
  user->scale1=4.08e0; user->scale2=0.417e0;

  info = PetscMalloc(sizeof(double)*user->mdiv4, (void*) &user->V0); CHKERRQ(info);
  info = PetscMalloc(sizeof(double)*user->mdiv4, (void*) &user->V1); CHKERRQ(info);
  //  user->V0 = (double*)PetscMalloc(sizeof(double)*user->mdiv4); CHKPTRQ(user->V0);
  //  user->V1 = (double*)PetscMalloc(sizeof(double)*user->mdiv4); CHKPTRQ(user->V1);

  info = PetscMalloc(sizeof(double)*user->mdiv2, (void*) &user->y); CHKERRQ(info);
  //  user->y = (double*)PetscMalloc(sizeof(double)*user->mdiv2); CHKPTRQ(user->y); 

  user->V0[0]=0.7140e0;  user->V0[1]=0.7169e0;  user->V0[2]=0.7232e0; 
  user->V0[3]=0.7151e0;  user->V0[4]=0.6848e0;  user->V0[5]=0.7070e0; 
  user->V0[6]=0.7177e0;  user->V0[7]=0.7073e0;  user->V0[8]=0.6734e0; 
  user->V0[9]=0.7174e0;  user->V0[10]=0.7125e0; user->V0[11]=0.6947e0; 
  user->V0[12]=0.7121e0; user->V0[13]=0.7166e0; user->V0[14]=0.6894e0;
  user->V0[15]=0.6897e0; user->V0[16]=0.7024e0; user->V0[17]=0.7026e0; 
  user->V0[18]=0.6800e0; user->V0[19]=0.6957e0; user->V0[20]=0.6987e0; 
  user->V0[21]=0.7111e0; user->V0[22]=0.7097e0; user->V0[23]=0.6809e0; 
  user->V0[24]=0.7139e0; user->V0[25]=0.7046e0; user->V0[26]=0.6950e0; 
  user->V0[27]=0.7032e0; user->V0[28]=0.7019e0; user->V0[29]=0.6975e0; 
  user->V0[30]=0.6955e0; user->V0[31]=0.7056e0; user->V0[32]=0.6965e0; 
  user->V0[33]=0.6848e0; user->V0[34]=0.6995e0; user->V0[35]=0.6105e0; 
  user->V0[36]=0.6027e0; user->V0[37]=0.6084e0; user->V0[38]=0.6081e0; 
  user->V0[39]=0.6057e0; user->V0[40]=0.6116e0; user->V0[41]=0.6052e0; 
  user->V0[42]=0.6136e0; user->V0[43]=0.6032e0; user->V0[44]=0.6081e0;
  user->V0[45]=0.6092e0; user->V0[46]=0.6122e0; user->V0[47]=0.6157e0; 
  user->V0[48]=0.6191e0; user->V0[49]=0.6169e0; user->V0[50]=0.5483e0; 
  user->V0[51]=0.5371e0; user->V0[52]=0.5576e0; user->V0[53]=0.5521e0; 
  user->V0[54]=0.5495e0; user->V0[55]=0.5499e0; user->V0[56]=0.4937e0; 
  user->V0[57]=0.5092e0; user->V0[58]=0.5433e0; user->V0[59]=0.5018e0; 
  user->V0[60]=0.5363e0; user->V0[61]=0.4977e0; user->V0[62]=0.5296e0;
  
  user->V1[0]=5.145e0;  user->V1[1]=5.241e0;  user->V1[2]=5.389e0; 
  user->V1[3]=5.211e0;  user->V1[4]=5.154e0;  user->V1[5]=5.105e0; 
  user->V1[6]=5.191e0;  user->V1[7]=5.013e0;  user->V1[8]=5.582e0; 
  user->V1[9]=5.208e0;  user->V1[10]=5.142e0; user->V1[11]=5.284e0; 
  user->V1[12]=5.262e0; user->V1[13]=6.838e0; user->V1[14]=6.215e0; 
  user->V1[15]=6.817e0; user->V1[16]=6.889e0; user->V1[17]=6.732e0; 
  user->V1[18]=6.717e0; user->V1[19]=6.468e0; user->V1[20]=6.776e0; 
  user->V1[21]=6.574e0; user->V1[22]=6.465e0; user->V1[23]=6.090e0; 
  user->V1[24]=6.350e0; user->V1[25]=4.255e0; user->V1[26]=4.154e0; 
  user->V1[27]=4.211e0; user->V1[28]=4.287e0; user->V1[29]=4.104e0; 
  user->V1[30]=4.007e0; user->V1[31]=4.261e0; user->V1[32]=4.150e0; 
  user->V1[33]=4.040e0; user->V1[34]=4.155e0; user->V1[35]=5.086e0; 
  user->V1[36]=5.021e0; user->V1[37]=5.040e0; user->V1[38]=5.247e0; 
  user->V1[39]=5.125e0; user->V1[40]=5.136e0; user->V1[41]=4.949e0; 
  user->V1[42]=5.253e0; user->V1[43]=5.154e0; user->V1[44]=5.227e0; 
  user->V1[45]=5.120e0; user->V1[46]=5.291e0; user->V1[47]=5.294e0; 
  user->V1[48]=5.304e0; user->V1[49]=5.209e0; user->V1[50]=5.384e0; 
  user->V1[51]=5.490e0; user->V1[52]=5.563e0; user->V1[53]=5.532e0; 
  user->V1[54]=5.372e0; user->V1[55]=5.423e0; user->V1[56]=7.237e0; 
  user->V1[57]=6.944e0; user->V1[58]=6.957e0; user->V1[59]=7.138e0; 
  user->V1[60]=7.009e0; user->V1[61]=7.074e0; user->V1[62]=7.046e0;


  user->y[0]=9.3636e0;    user->y[1]=9.3512e0;    user->y[2]=9.4891e0; 
  user->y[3]=9.1888e0;    user->y[4]=9.3161e0;    user->y[5]=9.2585e0;  
  user->y[6]=9.2913e0;    user->y[7]=9.3914e0;    user->y[8]=9.4524e0; 
  user->y[9]=9.4995e0;    user->y[10]=9.4179e0;   user->y[11]=9.468e0; 
  user->y[12]=9.4799e0;   user->y[13]=11.2917e0;  user->y[14]=11.5062e0; 
  user->y[15]=11.4579e0;  user->y[16]=11.3977e0;  user->y[17]=11.3688e0;
  user->y[18]=11.3897e0;  user->y[19]=11.3104e0;  user->y[20]=11.3882e0;
  user->y[21]=11.3629e0;  user->y[22]=11.3149e0;  user->y[23]=11.2474e0; 
  user->y[24]=11.2507e0;  user->y[25]=8.1678e0;   user->y[26]=8.1017e0;  
  user->y[27]=8.3506e0;   user->y[28]=8.3651e0;   user->y[29]=8.2994e0; 
  user->y[30]=8.1514e0;   user->y[31]=8.2229e0;   user->y[32]=8.1027e0; 
  user->y[33]=8.3785e0;   user->y[34]=8.4118e0;   user->y[35]=8.0955e0; 
  user->y[36]=8.0613e0;   user->y[37]=8.0979e0;   user->y[38]=8.1364e0;  
  user->y[39]=8.1700e0;   user->y[40]=8.1684e0;   user->y[41]=8.0885e0; 
  user->y[42]=8.1839e0;   user->y[43]=8.1478e0;   user->y[44]=8.1827e0;  
  user->y[45]=8.029e0;    user->y[46]=8.1000e0;   user->y[47]=8.2579e0; 
  user->y[48]=8.2248e0;   user->y[49]=8.2540e0;   user->y[50]=6.8518e0;  
  user->y[51]=6.8547e0;   user->y[52]=6.8831e0;   user->y[53]=6.9137e0; 
  user->y[54]=6.8984e0;   user->y[55]=6.8888e0;   user->y[56]=8.5189e0; 
  user->y[57]=8.5308e0;   user->y[58]=8.5184e0;   user->y[59]=8.5222e0; 
  user->y[60]=8.5705e0;   user->y[61]=8.5353e0;   user->y[62]=8.5213e0;  
  user->y[63]=8.3158e0;   user->y[64]=8.1995e0;   user->y[65]=8.2283e0;
  user->y[66]=8.1857e0;   user->y[67]=8.2738e0;   user->y[68]=8.2131e0; 
  user->y[69]=8.2613e0;   user->y[70]=8.2315e0;   user->y[71]=8.2078e0; 
  user->y[72]=8.2996e0;   user->y[73]=8.3026e0;   user->y[74]=8.0995e0;  
  user->y[75]=8.2990e0;   user->y[76]=9.6753e0;   user->y[77]=9.6687e0; 
  user->y[78]=9.5704e0;   user->y[79]=9.5435e0;   user->y[80]=9.6780e0;
  user->y[81]=9.7668e0;   user->y[82]=9.7827e0;   user->y[83]=9.7844e0; 
  user->y[84]=9.7011e0;   user->y[85]=9.8006e0;   user->y[86]=9.7610e0;  
  user->y[87]=9.7813e0;   user->y[88]=7.3073e0;   user->y[89]=7.2572e0; 
  user->y[90]=7.4686e0;   user->y[91]=7.3659e0;   user->y[92]=7.3587e0;  
  user->y[93]=7.3132e0;   user->y[94]=7.3542e0;   user->y[95]=7.2339e0; 
  user->y[96]=7.4375e0;   user->y[97]=7.4022e0;   user->y[98]=10.7914e0; 
  user->y[99]=10.6554e0;  user->y[100]=10.7359e0; user->y[101]=10.7583e0;
  user->y[102]=10.7735e0; user->y[103]=10.7907e0; user->y[104]=10.6465e0;
  user->y[105]=10.6994e0; user->y[106]=10.7756e0; user->y[107]=10.7402e0;
  user->y[108]=10.6800e0; user->y[109]=10.7000e0; user->y[110]=10.8160e0;
  user->y[111]=10.6921e0; user->y[112]=10.8677e0; user->y[113]=12.3495e0;
  user->y[114]=12.4424e0; user->y[115]=12.4303e0; user->y[116]=12.5086e0;
  user->y[117]=12.4513e0; user->y[118]=12.4625e0; user->y[119]=16.2290e0;
  user->y[120]=16.2781e0; user->y[121]=16.2082e0; user->y[122]=16.2715e0;
  user->y[123]=16.2464e0; user->y[124]=16.1626e0; user->y[125]=16.1568e0;


  return 0;
}


int AppCtxDestroy(AppCtx *user)
{

  PetscFree(user->V0);
  PetscFree(user->V1);
  PetscFree(user->y);

  return 0;
}


