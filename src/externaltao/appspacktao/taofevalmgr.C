
#include "fevalmgr.H"
#include "gci.H"
#include "tao_solver.h"
#include "taofevalwkr.h"
#include "taofevalmgr.h"


// THIS SPECIAL FUNCTINO MUST BE DEFINED AFTER CLASS MYFEVALMGR. The
// APPS object files never know the actual MyFevalMgr object, only its
// base class FevalMgr. A new myfevalmgr.o object file can be linked
// in without recompiling anything else. The secret is this special
// funcion which returns a pointer to the base class.  -- From APPS


#undef __FUNCT__
#define __FUNCT__ "myfevalmgr"
// *** This MUST be called 'myfevalmgr' ***
FevalMgr* myfevalmgr()
{
  return new TaoFevalMgr();
}

//==========================================
TAO_SOLVER TaoFevalMgr::tao = 0;
TAO_APPS* TaoFevalMgr::appsPtr = 0;


//==========================================
#undef __FUNCT__
#define __FUNCT__ "TaoFevalMgr::TaoFevalMgr"
TaoFevalMgr::TaoFevalMgr()
{
  this->xl=0;
  this->xu=0;
}


//==========================================
#undef __FUNCT__
#define __FUNCT__ "TaoFevalMgr::~TaoFevalMgr"
TaoFevalMgr::~TaoFevalMgr() {
  TaoVecDestroy(this->xl);
  TaoVecDestroy(this->xu);
}


//==========================================
#undef __FUNCT__
#define __FUNCT__ "TaoFevalMgr::setTao"
int TaoFevalMgr::setTao(TAO_SOLVER t, TAO_APPS *ap)
{
  int info;
  TaoFunctionBegin;

  tao = t;
  appsPtr = ap;
  TaoFunctionReturn(0);
}

//==========================================
#undef __FUNCT__
#define __FUNCT__ "~TaoFevalMgr"
bool TaoFevalMgr::setFevalMgr(const vector<string>& fargs, bool issingle)
{
  int info;

  // This sets up all the parameters for the function evaluation using
  // the apps command line arguments 2 and up. Here 'argc' = number of
  // fevalMgr params and 'argv' = array of parameters.

  // tao must be set before calling this function
  if (!tao)
  {
    // Should find a way to give a reason for exiting
    return false;
  }

  // Get the problem size and should check that it is valid
  ndim = appsPtr->ndim;

  
  setBounds();

  // The feval is the manager of the parallel operations if issingle
  // is true.
  ismanager = issingle;

  // Set initx 
  setInitX();


  return true;
} 

//==========================================
#undef __FUNCT__
#define __FUNCT__ "TaoFevalMgr::setbounds"
void TaoFevalMgr::setBounds(void)
{
  TaoVec *LL, *UU;
  double *l, *u;
  int i, dim;
  int info;

  
  // Get the bounds from tao
  info = TaoGetVariableBounds(tao, &LL, &UU);
  if (LL && UU)
  {
    usebounds = true;
    info = LL->GetArray(&l,&dim);
    info = UU->GetArray(&u,&dim);
  
    // Convert l, u, to STP vectors lower, upper
    lower.reserve(dim);
    upper.reserve(dim);
    for (i=0;i<dim;i++)
    {
      lower.push_back(l[i]);
      upper.push_back(u[i]);
    }

    info = LL->RestoreArray(&l,&dim);
    info = UU->RestoreArray(&u,&dim);

  }
  else
    usebounds = false;

  
}

//==========================================
#undef __FUNCT__
#define __FUNCT__ "TaoFevalMgr::setInitX"
void TaoFevalMgr::setInitX(void)
{
  TaoVec *XX0;
  double *x0;
  int i;
  int info,dim;

  // Get the initial value of the vector from tao
  info = TaoGetSolution(tao,&XX0);
  info = XX0->GetArray(&x0,&dim);

  initx.resize(ndim);
  // the initial vector must be stored in the class member STL vector initx 
  for (i=0;i<ndim;i++)
    initx[i] = x0[i];

  info = XX0->RestoreArray(&x0,&dim);
  
}

