/*$Id$*/

/* Program usage: mpirun -np 1 rosenbrock1 [-help] [all TAO options] */

/* T 
   Concepts: TAO - Solving an unconstrained minimization problem
   Routines: TaoABCApplication::SetDimension(); TaoABCApplication::GetSolutionAndGradient();
   Routines: TaoABCApplication::ComputeObjectiveAndGradient(); 
   Processors: 1
T*/ 

//   Routines: TaoABCApplication::InitializeVariables(); TaoABCApplication::ComputeVariableBounds();


static  char help[] = 
"This example demonstrates use of the TAO package to \n\
solve an unconstrained minimization problem on a single processor.  We \n\
minimize the extended Rosenbrock function: \n\
   sum_{i=0}^{n/2-1} ( alpha*(x_{2i+1}-x_{2i}^2)^2 + (1-x_{2i})^2 ) \n";

/* ---------------------------------------------------------------------- */

#include "tao_solver.h"
#include <math.h>

class Rosenbrock: public TaoABCApplication{
 private:
  double alpha;
 public:
  Rosenbrock();
  int ComputeObjectiveAndGradient(double*, int, double*, double *);
};

#undef __FUNCT__  
#define __FUNCT__ "Rosenbrock::Rosenbrock"
Rosenbrock::Rosenbrock(){
  this->alpha=99;
  this->SetNumberOfVariables(2);
  return;
}

#undef __FUNCT__
#define __FUNCT__ "Rosenbrock::ComputeObjectiveAndGradient"
int Rosenbrock::ComputeObjectiveAndGradient(double *x, int n, double *f, double *g)
{
  int    i,nn=n/2;
  double t1,t2,ff=0;

  for (i=0; i<nn; i++){
    t1 = x[2*i+1]-x[2*i]*x[2*i]; t2= 1-x[2*i];
    ff += alpha*t1*t1 + t2*t2;
    g[2*i] = -4*alpha*t1*x[2*i]-2.0*t2;
    g[2*i+1] = 2*alpha*t1;
  }
  *f=ff;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  TAO_SOLVER tao;                   /* TAO_SOLVER solver context */
  Rosenbrock *rose;
  double     *x;
  int        i,n, info;                  /* error code */

  /* Initialize TAO */
  //  info=PetscInitialize(&argc,&argv,(char *)0,help); CHKERRQ(info);
  info=TaoInitialize(&argc,&argv,(char *)0,help); CHKERRQ(info);

  info = TaoCreate(MPI_COMM_SELF,"tao_lmvm",&tao); CHKERRQ(info);
  rose = new Rosenbrock();
  
  info = TaoSetApplication(tao,rose); CHKERRQ(info);
  info = TaoSetFromOptions(tao);CHKERRQ(info);

  info = TaoSolve(tao); CHKERRQ(info);

  /* Print Solution */
  info = rose->GetSolution(x, n); CHKERRQ(info);
  printf("\nSolution:  ");
  for (i=0; i<n; i++){ printf(" %7.4e ",x[i]); }
  printf("\n");

  /* Free data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoDestroyApplication(rose); CHKERRQ(info);

  /* Finalize TAO */
  info = TaoFinalize(); CHKERRQ(info);
  //  info = PetscFinalize(); CHKERRQ(info);
  return 0;
}

