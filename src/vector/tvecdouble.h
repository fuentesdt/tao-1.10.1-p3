#ifndef SIMPLEVECTOR_H
#define SIMPLEVECTOR_H

// #include "tao_basictypes.h"
#include "taovec.h"

class TaoVecDoubleArray: public TaoVec{

 protected:

  TaoInt n;
  double *v;
  TaoInt dallocated;
 public:
  
  inline int GetData(double**dd,TaoInt*nn){*dd=v;*nn=n; return 0;}

  TaoVecDoubleArray( TaoInt nn );
  TaoVecDoubleArray( TaoInt nn , double *vv);
  ~TaoVecDoubleArray(){if (n>0 && dallocated) delete [] v;};


  int Compatible (TaoVec *v, TaoTruth*);
  int GetArray(TaoScalar **, TaoInt*);
  int RestoreArray(TaoScalar **, TaoInt*);
  int Clone(TaoVec**);
  int GetDimension(TaoInt *);

  int View();

  int GetDoubles(double **, TaoInt*);
  int RestoreDoubles(double **, TaoInt*);

};

#endif

