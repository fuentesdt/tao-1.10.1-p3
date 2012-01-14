#include "tao_general.h"
#include "taois.h"

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetDestroy"
/*@C
   TaoIndexSetDestroy - Destroys the TaoIndexSet object.

   Input Parameter:
.  SS - the vector

   Level: beginner
@*/
int TaoIndexSetDestroy( TaoIndexSet * SS){
  TaoFunctionBegin;
  if (SS) delete SS;
  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::Duplicate"
/*@C
   Duplicate - Creates a new TaoIndexSet object with the same structure as this one.

   Output Parameter:
.  SS -  new TaoIndexSet object

   Level: intermediate
@*/
int TaoIndexSet::Duplicate(TaoIndexSet** SS){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}


#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::IsSame"
/*@C
   IsSame - Determines whether this index set is the same as another.

   Input Parameter:
.  ss -  an index set

   Output Parameter:
.  flg -  set to TAO_TRUE if the two index sets are the same and TAO_FALSE otherwise

   Level: intermediate
@*/
int TaoIndexSet::IsSame(TaoIndexSet* ss, TaoTruth* flg){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::ComplementOf"
/*@C
   ComplementOf - Let this index set be the complement of the input index set.

   Input Parameter:
.  ss -  an index set

   Level: intermediate
@*/
int TaoIndexSet::ComplementOf(TaoIndexSet* ss){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::IntersectionOf"
/*@C
   IntersectionOf - Let this index set be the intersection of the given two index sets.

   Input Parameter:
+  ss1 -  an index set
-  ss2 -  an index set

   Level: intermediate
@*/
int TaoIndexSet::IntersectionOf(TaoIndexSet* ss1, TaoIndexSet* ss2){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::UnionOf"
/*@C
   UnionOf - Let this index set be the union of the given two index sets.

   Input Parameter:
+  ss1 -  an index set
-  ss2 -  an index set

   Level: intermediate
@*/
int TaoIndexSet::UnionOf(TaoIndexSet* ss1, TaoIndexSet* ss2){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::WhichEqual"
/*@C
   WhichEqual - Describes which elements of the vectors equal one another.

   Input Parameter:
.  vv1, vv2 -  the vectors

   Level: intermediate
@*/
int TaoIndexSet::WhichEqual(TaoVec* vv1,TaoVec* vv2){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

int TaoIndexSet::WhichLessThan(TaoVec* tv1,TaoVec* tv2){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

int TaoIndexSet::WhichGreaterThan(TaoVec* tv1,TaoVec* tv2){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::WhichBetween"
/*@C
   WhichBetween - Describes which elements of vector tv are greater than
   the corresponding element of vector low and less than corresponding element
   of vector high.

   Input Parameters:
+  xxll -  the vector representing lower bounds on the vector tv
.  vv - the TAO vector
-  xxuu - the vector of upper bounds

   Level: intermediate
@*/
int TaoIndexSet::WhichBetween(TaoVec* xxll,TaoVec* vv ,TaoVec* xxuu){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::WhichBetweenOrEqual"
/*@C
   WhichBetweenOrEqual - Describes which elements of vector are in the closed
   interface defined by the lower and upper bounds.

   Input Parameters:
+  xxll -  the vector representing lower bounds on the vector tv
.  vv - the TAO vector
-  xxuu - the vector of upper bounds

   Level: intermediate
@*/
int TaoIndexSet::WhichBetweenOrEqual(TaoVec *xxll, TaoVec *vv, TaoVec *xxuu)
{
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::GetSize"
/*@C
   GetSize - Describes which elements of vector tv are greater than
   the corresponding element of vector low and less than corresponding element
   of vector high.

   Output Parameter:
.  nn - the number of size of the index set.

   Level: intermediate
@*/
int TaoIndexSet::GetSize(TaoInt *nn){;
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}


#undef __FUNCT__
#define __FUNCT__ "TaoIndexSet::View"
/*@C
  View - Views the contents of the index set.
  
  Input Parameters: none
  
  Level: intermediate
@*/
int TaoIndexSet::View(){
  TaoFunctionBegin;
  SETERRQ(56,"Operation not defined");
  /* TaoFunctionReturn(1); */
}

