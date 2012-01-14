#include "tao_general.h"            /*I "tao.h" I*/
#ifdef TAO_USE_PETSC
#include "taois_petsc.h"
#include "../vector/taovec_petsc.h"

int VecWhichBetween(Vec, Vec, Vec, IS *);
int VecWhichBetweenOrEqual(Vec, Vec, Vec, IS *);
int VecWhichGreaterThan(Vec, Vec, IS * );
int VecWhichLessThan(Vec, Vec, IS *);
int VecWhichEqual(Vec, Vec, IS *);

static PetscLogEvent petscisevents=0;
static PetscLogEvent tredistribute=0;
//static TaoPetscISType defaultistype=TaoRedistributeSubset;
//static TaoPetscISType defaultistype=TaoNoRedistributeSubset;
static TaoPetscISType defaultistype=TaoMaskFullSpace;

#undef __FUNCT__
#define __FUNCT__ "TaoWrapPetscIS"
/* @C
   TaoWrapPetscIS - Create a new TaoIndexSet object using PETSc.

   Input Parameter:
+  S -  a PETSc IS
-  imax - the maximum local length of the index set (Should be equal to or greater than the local length of the vector)

   Output Parameter:
.  SS - address of a new TaoIndexSet

   Note:  
   A TaoIndexSetPetsc is an object with the methods of an abstract
   TaoIndexSet object.  A TaoIndexSetPetsc contains an implementation of the 
   TaoIndexSet methods.  Routines using these vectors should declare a 
   pointer to a TaoIndexSet, assign this pointer to the address of a 
   TaoIndexSet object,use the pointer to invoke methods on the object, and use
   this pointer as an argument when calling other routines.  This usage is 
   different from the usage of a PETSc IS.  In PETSc, applications will 
   typically declare a IS, and pass it as an argument into routines.  That is,
   applications will typically declare a pointer to a TaoIndexSet and use the
   pointer, or declare a IS and use it directly.

.seealso TaoIndexSetGetPetscIS(), TaoIndexSetDestroy()

   Level: developer
@ */
int TaoWrapPetscIS( IS S, PetscInt imax, TaoIndexSetPetsc ** SS){
  PetscFunctionBegin;
  PetscValidHeaderSpecific(S,IS_COOKIE,1);
  //  if (0==0) return 1;
  //  *SS = new  TaoIndexSetPetsc(imax, S);
  PetscFunctionReturn(1);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::TaoIndexSetPetsc"
/* @C
   TaoIndexSetPetsc - Create a new TaoIndexSet object using PETSc.

   Input Parameter:
-  V = A vector from the space that this indexset will describe.
+  SS - an Index Set

   Level: beginner

   Note: 
   This index set will be deleted when the TaoIndexSet is deleted

   Note: 
   The constuctor can also be called without the IndexSet
   
   Note:
   The method TaoIndexSetPetsc::GetIS() returns a pointer to the current PETSc Index Set.
@ */
TaoIndexSetPetsc::TaoIndexSetPetsc(Vec V, IS SS):TaoIndexSet(){
  Vec V2;
  isp=0;isp2=0;ISGathered=0;ispComplement=0;scatter=0;
  VecDuplicate(V,&V2);
  this->SetUp(V2,SS);
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::TaoIndexSetPetsc"
/* @C
   TaoIndexSetPetsc - Create a new TaoIndexSet object using PETSc.

   Input Parameter:
.  V = A vector from the space that this indexset will describe.

   Level: beginner

   Note:
   The method TaoIndexSetPetsc::GetIS() returns a pointer to the current PETSc Index Set.
@ */
TaoIndexSetPetsc::TaoIndexSetPetsc(Vec V):TaoIndexSet(){
  int info;
  PetscInt low,high;
  IS is;
  Vec V2;
  isp=0;isp2=0;ISGathered=0;ispComplement=0;scatter=0;
  VecDuplicate(V,&V2);
  info=VecGetOwnershipRange(V,&low,&high);
  info=ISCreateStride(PETSC_COMM_WORLD,1,low,1,&is);
  this->SetUp(V2,is);
  info = PetscObjectDereference((PetscObject)is);
  return;
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::~TaoIndexSetPetsc"
TaoIndexSetPetsc::~TaoIndexSetPetsc(){
  clearit();
  delete[] iptr; 
  if (this->VSpace){
    VecDestroy(this->VSpace);
  }
};

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::SetUp"
int TaoIndexSetPetsc::SetUp(Vec V2, IS SS){

  int info;
  int size;
  PetscInt imax;
  PetscTruth flg,ivalue;
  MPI_Comm comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V2,VEC_COOKIE,1);

  info=VecGetLocalSize(V2,&imax); CHKERRQ(info);
  VSpace=V2;

  iptr = new PetscInt[imax];

  isp=0;isp2=0;ISGathered=0;ispComplement=0;scatter=0;
  this->ispviewer=PETSC_VIEWER_STDOUT_WORLD;
  reducedtype=defaultistype;

  info=this->SetIS(SS); CHKERRQ(info);
  nlocal=imax;

  info = PetscObjectGetComm((PetscObject)V2,&comm);CHKERRQ(info);
  MPI_Comm_size(comm,&size);
  if (size==1 && defaultistype!=TaoMaskFullSpace &&defaultistype!=TaoMatrixFree ){ 
    reducedtype=TaoSingleProcessor;
  } else {
    PetscOptionsGetTruth(0,"-redistribute",&ivalue,&flg);
    if (flg && ivalue==PETSC_FALSE) reducedtype=TaoNoRedistributeSubset;
    else if (flg && ivalue==PETSC_TRUE)reducedtype=TaoRedistributeSubset;
    else reducedtype=defaultistype;
  }
  PetscOptionsHasName(0,"-taomask",&flg);
  if (flg) { reducedtype= TaoMaskFullSpace; }
  PetscOptionsHasName(0,"-taosubmatrixfree",&flg);
  if (flg) { reducedtype= TaoMatrixFree; }

  if (petscisevents==0){
    PetscLogEventRegister("Identify Indices",IS_COOKIE,&petscisevents);
  }
  if (tredistribute==0){
    PetscLogEventRegister("IS Redistribute ",IS_COOKIE,&tredistribute);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::SetIS"
int TaoIndexSetPetsc::SetIS(IS SS){

  int info;
  PetscFunctionBegin;
  if (SS){
    PetscValidHeaderSpecific(SS,IS_COOKIE,1);
    PetscObjectReference((PetscObject)SS); 
  }
  info=this->clearit(); CHKERRQ(info);
  this->isp=SS;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoSelectSubset"
/*@C
   TaoSelectSubset - Select the manner in which subsets of variables in a matrix are selected

   Input Parameter:
.  type - the manner is which subsets are chosen.  

   Choice of types:
+  TaoRedistributeSubset - When using multiple processors, redistribute the subset of variables
over the processors.
.  TaoNoRedistributeSubset - When using multiple processors, keep variables on current processors
.  TaoMaskFullSpace - Use the full set of variables, and mask the unneeded ones.
-  TaoMatrixFree - Use the full set of variables, and mask the unneeded ones.

   Note:  
   For use with active set methods in bound constrained optimization such as TRON and GPCG

   Level: advanced
@*/
int TaoSelectSubset( TaoPetscISType type){
  PetscFunctionBegin;
  defaultistype=type;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::GetMask"
int TaoIndexSetPetsc::GetMask(Vec *VMask){
  int info;
  PetscScalar zero=0.0;

  PetscFunctionBegin;
  info=VecSet(this->VSpace, zero); CHKERRQ(info);
  info=VecISSetToConstant(this->isp,1.0,this->VSpace); CHKERRQ(info);
  *VMask=this->VSpace;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::UnionOf"
int TaoIndexSetPetsc::UnionOf(TaoIndexSet* S1, TaoIndexSet* S2){
  IS SS1=((TaoIndexSetPetsc*)S1)->GetIS();
  IS SS2=((TaoIndexSetPetsc*)S2)->GetIS();
  int info;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  info=ISSum(SS1,SS2,&SS1); CHKERRQ(info);
  info=PetscLogEventEnd(petscisevents,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::IntersectionOf"
int TaoIndexSetPetsc::IntersectionOf(TaoIndexSet* S1, TaoIndexSet* S2){
  IS SS1=((TaoIndexSetPetsc*)S1)->GetIS();
  IS SS2=((TaoIndexSetPetsc*)S2)->GetIS();
  IS SS3;
  int info;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  info=ISDifference(SS1,SS2,&SS3); CHKERRQ(info);
  info=this->SetIS(SS3); CHKERRQ(info);
  info=ISDestroy(SS3); CHKERRQ(info);
  info=PetscLogEventEnd(petscisevents,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::ComplementOf"
int TaoIndexSetPetsc::ComplementOf(TaoIndexSet* SS){
  IS SS1=((TaoIndexSetPetsc*)SS)->GetIS();
  IS SS2;
  int info;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  info=ISCreateComplement(SS1,this->VSpace,&SS2); CHKERRQ(info);
  info=this->SetIS(SS2); CHKERRQ(info);
  info=ISDestroy(SS2); CHKERRQ(info);
  info=PetscLogEventEnd(petscisevents,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::Duplicate"
int TaoIndexSetPetsc::Duplicate(TaoIndexSet**SS){
  int info;
  TaoIndexSetPetsc *S;
  IS is;

  PetscFunctionBegin;
  if (isp){
    info = ISDuplicate(isp,&is); CHKERRQ(info);
    S = new TaoIndexSetPetsc(this->VSpace, is);
    info=ISDestroy(is); CHKERRQ(info);
  } else {
    S = new TaoIndexSetPetsc(this->VSpace);
  }
  *SS=S;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::IsSame"
int TaoIndexSetPetsc::IsSame(TaoIndexSet*SS, TaoTruth*flg){
  int info;
  PetscTruth pflg;
  IS SS2=((TaoIndexSetPetsc*)SS)->GetIS();
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  info=ISEqual(isp,SS2,&pflg); CHKERRQ(info);
  if (pflg==PETSC_TRUE) *flg=TAO_TRUE; else *flg=TAO_FALSE;
  info=PetscLogEventEnd(petscisevents,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::View"
int TaoIndexSetPetsc::View(){
  int info;
  PetscFunctionBegin;
  info=ISView(this->isp,this->ispviewer);CHKERRQ(info); CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::WhichEqual"
int TaoIndexSetPetsc::WhichEqual(TaoVec*v1,TaoVec*v2){
  int info;
  TaoVecPetsc *tv1, *tv2;
  IS ISnew;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  tv1 = dynamic_cast<TaoVecPetsc*>(v1);
  tv2 = dynamic_cast<TaoVecPetsc*>(v2);
  if (tv1 == NULL) {
    SETERRQ(1,"Could not cast argument 1 to TaoVecPetsc*");
  }
  if (tv2 == NULL) {
    SETERRQ(1,"Could not cast argument 2 to TaoVecPetsc*");
  }
  info=VecWhichEqual(tv1->GetVec(),tv2->GetVec(),&ISnew); CHKERRQ(info);
  info=this->SetIS(ISnew); CHKERRQ(info);
  info=ISDestroy(ISnew); CHKERRQ(info);
  info=PetscLogEventEnd(petscisevents,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::WhichLessThan"
int TaoIndexSetPetsc::WhichLessThan(TaoVec*v1,TaoVec*v2){
  int info;
  Vec V1=((TaoVecPetsc *)v1)->GetVec();
  Vec V2=((TaoVecPetsc *)v2)->GetVec();
  IS ISnew;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  info=VecWhichLessThan(V1,V2,&ISnew); CHKERRQ(info);
  info=this->SetIS(ISnew); CHKERRQ(info);
  info=ISDestroy(ISnew); CHKERRQ(info);
  info=PetscLogEventEnd(petscisevents,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::WhichGreaterThan"
int TaoIndexSetPetsc::WhichGreaterThan(TaoVec*v1,TaoVec*v2){
  int info;
  Vec V1=((TaoVecPetsc *)v1)->GetVec();
  Vec V2=((TaoVecPetsc *)v2)->GetVec();
  IS ISnew;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  info=VecWhichGreaterThan(V1,V2,&ISnew); CHKERRQ(info);
  info=this->SetIS(ISnew); CHKERRQ(info);
  info=ISDestroy(ISnew); CHKERRQ(info);
  info=PetscLogEventEnd(petscisevents,0,0,0,0); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::WhichBetween"
int TaoIndexSetPetsc::WhichBetween(TaoVec*V1,TaoVec*V2,TaoVec*V3){
  int info;
  Vec VecLow=((TaoVecPetsc *)V1)->GetVec();
  Vec VecHigh=((TaoVecPetsc *)V3)->GetVec();
  Vec V=((TaoVecPetsc *)V2)->GetVec();
  IS ISnew;

  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  info=VecWhichBetween(VecLow,V,VecHigh,&ISnew); CHKERRQ(info);
  info=this->SetIS(ISnew); CHKERRQ(info);
  info=ISDestroy(ISnew); CHKERRQ(info);
  info=PetscLogEventEnd(petscisevents,0,0,0,0);CHKERRQ(info);          

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::WhichBetweenOrEqual"
int TaoIndexSetPetsc::WhichBetweenOrEqual(TaoVec *V1, TaoVec *V2, TaoVec *V3)
{
  int info;
  Vec VecLow  = ((TaoVecPetsc *)V1)->GetVec();
  Vec VecHigh = ((TaoVecPetsc *)V3)->GetVec();
  Vec V       = ((TaoVecPetsc *)V2)->GetVec();
  IS ISnew;

  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0); CHKERRQ(info);
  info=VecWhichBetweenOrEqual(VecLow,V,VecHigh,&ISnew); CHKERRQ(info);
  info=this->SetIS(ISnew); CHKERRQ(info);
  info=ISDestroy(ISnew); CHKERRQ(info);
  info=PetscLogEventEnd(petscisevents,0,0,0,0);CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::GetSize"
int TaoIndexSetPetsc::GetSize(TaoInt *nn){
  PetscFunctionBegin;
  int info=ISGetSize(isp,nn); CHKERRQ(info);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::clearit"
int TaoIndexSetPetsc::clearit(){
  int info;
  PetscFunctionBegin;
  if (scatter){
    info=VecScatterDestroy(scatter); CHKERRQ(info);
    scatter=0;
  }
  if (ISGathered){
    info=ISDestroy(ISGathered); CHKERRQ(info);
    ISGathered=0;
  }
  if (isp){
    info=ISDestroy(isp); CHKERRQ(info);
    isp=0;
  }
  if (isp2){
    info=ISDestroy(isp2); CHKERRQ(info);
    isp2=0;
  }
  if (ispComplement){
    info=ISDestroy(ispComplement); CHKERRQ(info);
    ispComplement=0;
  }
  isp=0;isp2=0;ISGathered=0;ispComplement=0;scatter=0;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::Redistribute"
int TaoIndexSetPetsc::RedistributeIS(IS*ISNew){
  int info;
  PetscInt n,nn,low,high,i;
  const PetscInt *is;
  PetscInt *gl;
  Vec R;
  IS ISAll;
  MPI_Comm comm;

  PetscFunctionBegin;

  info=PetscLogEventBegin(petscisevents,0,0,0,0);CHKERRQ(info);
  info = VecGetSize(VSpace,&n); CHKERRQ(info);
  info = this->GetSize(&nn); CHKERRQ(info);
  if (reducedtype==TaoMaskFullSpace){
    *ISNew = isp;
  } else if (reducedtype==TaoSingleProcessor){
    *ISNew = isp;
  } else if (reducedtype==TaoNoRedistributeSubset || n==nn){
    *ISNew = isp;
  } else if (reducedtype==TaoRedistributeSubset){

    info=PetscLogEventBegin(tredistribute,0,0,0,0);CHKERRQ(info);

    if (isp2==0){
      info = this->GetWholeIS(&ISAll); CHKERRQ(info);
      
      info = PetscObjectGetComm((PetscObject)isp,&comm);CHKERRQ(info);
      info = VecCreateMPI(comm,PETSC_DECIDE,nn,&R);CHKERRQ(info);
      info = VecGetLocalSize(R,&n); CHKERRQ(info);
      info = VecGetOwnershipRange(R, &low, &high); CHKERRQ(info);
      info = ISGetIndices(ISAll, &is); CHKERRQ(info);
      if (n>0){
	info = PetscMalloc( n*sizeof(PetscInt),&gl ); CHKERRQ(info);
	for (i=0; i<n; i++){
	  gl[i]= is[low+i];
	} 
      } else gl=NULL;
      
      info = PetscObjectGetComm((PetscObject)isp,&comm);CHKERRQ(info);
      info = ISCreateGeneral(comm,n,gl,ISNew); CHKERRQ(info);
      info = ISRestoreIndices(ISGathered, &is); CHKERRQ(info);
      if (n>0) {
	info = PetscFree(gl); CHKERRQ(info);
      }

      isp2=*ISNew;
      info = VecDestroy(R); CHKERRQ(info);

    }
    info=PetscLogEventEnd(tredistribute,0,0,0,0);CHKERRQ(info);
    *ISNew = isp2;
    
  }

  info=PetscLogEventEnd(petscisevents,0,0,0,0);CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::GetWholeIS"
int TaoIndexSetPetsc::GetWholeIS(IS*ISAll){
  int info;
  PetscFunctionBegin;
  info=PetscLogEventBegin(petscisevents,0,0,0,0);CHKERRQ(info);
  if (reducedtype==TaoSingleProcessor){
    *ISAll = isp;
  } else {
    info=PetscLogEventBegin(tredistribute,0,0,0,0);CHKERRQ(info);
    if (ISGathered==0){
      info = ISAllGather(isp,&ISGathered); CHKERRQ(info);
      info = ISSort(ISGathered); CHKERRQ(info);
    } else if (ISGathered==0){
      PetscInt  N;
      PetscMPIInt *sizes,*offsets,size,mpin;
      const PetscInt *lindices;
      PetscInt *indices;
      PetscInt i,n;
      MPI_Comm comm;

      info=PetscObjectGetComm((PetscObject)isp,&comm);CHKERRQ(info);
      MPI_Comm_size(comm,&size);
      PetscMalloc(2*size*sizeof(PetscInt),&sizes);
      offsets = sizes + size;

      info=ISGetLocalSize(isp,&n);CHKERRQ(info);
      mpin = (PetscMPIInt) n;
      if ((PetscInt)mpin != n) {
	SETERRQ(1,"Did not handle 64-bit-index correctly");
      }
      MPI_Allgather(&mpin,1,MPI_INT,sizes,1,MPI_INT,comm);
      offsets[0] = 0;
      for (i=1;i<size; i++) offsets[i] = offsets[i-1] + sizes[i-1];
      N = offsets[size-1] + sizes[size-1];

      PetscMalloc((N+1)*sizeof(PetscInt),&indices);
      info=ISGetIndices(isp,&lindices);CHKERRQ(info);
      MPI_Allgatherv((void*)lindices,n,MPIU_INT,indices,sizes,offsets,MPIU_INT,comm);
      info=ISRestoreIndices(isp,&lindices);CHKERRQ(info);

      info=ISCreateGeneral(comm,N,indices,&ISGathered);CHKERRQ(info);
      info = PetscFree(indices); CHKERRQ(info);
      info = PetscFree(sizes);  CHKERRQ(info);
    }
    *ISAll = ISGathered;
    info=PetscLogEventEnd(tredistribute,0,0,0,0);CHKERRQ(info);
  }

  info=PetscLogEventEnd(petscisevents,0,0,0,0);CHKERRQ(info);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::GetReducedVecScatter"
int TaoIndexSetPetsc::GetReducedVecScatter(Vec VFull, Vec VSub, VecScatter*scatterit){
  int info;
  PetscInt n,nlow,nhigh;
  IS S1, Ident;
  MPI_Comm comm;

  PetscFunctionBegin;

  if (scatter==0){
    info = VecGetSize(VFull,&n); CHKERRQ(info);
    info = VecGetSize(VSub,&nlow); CHKERRQ(info);
    if (reducedtype==TaoMaskFullSpace || n==nlow){

      info = VecScatterCreate(VFull,isp,VSub,isp,&scatter); CHKERRQ(info);

    } else if (reducedtype==TaoSingleProcessor || reducedtype==TaoNoRedistributeSubset){

      S1=isp;
      info = VecGetOwnershipRange(VSub,&nlow,&nhigh); CHKERRQ(info);
      n=nhigh-nlow;

      info = PetscObjectGetComm((PetscObject)VFull,&comm);CHKERRQ(info);
      info = ISCreateStride(comm,n,nlow,1,&Ident); CHKERRQ(info);
      info = VecScatterCreate(VFull,S1,VSub,Ident,&scatter); CHKERRQ(info);
      info = ISDestroy(Ident); CHKERRQ(info);

    } else {

      info = this->GetWholeIS(&S1); CHKERRQ(info);
      nlow = 0;
      info = VecGetSize(VSub,&n);CHKERRQ(info);

      info = PetscObjectGetComm((PetscObject)VFull,&comm);CHKERRQ(info);
      info = ISCreateStride(comm,n,nlow,1,&Ident); CHKERRQ(info);  
      info = VecScatterCreate(VFull,S1,VSub,Ident,&scatter); CHKERRQ(info);
      info = ISDestroy(Ident); CHKERRQ(info);

    }
  }
  
  *scatterit=scatter;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::GetReducedType"
int TaoIndexSetPetsc::GetReducedType(TaoPetscISType* rtype){
  PetscFunctionBegin;
  *rtype=this->reducedtype;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TaoIndexSetPetsc::GetComplementIS"
int TaoIndexSetPetsc::GetComplementIS(IS*isppp){
  int info;
  PetscFunctionBegin;
  if (this->ispComplement==0){
    info=PetscLogEventBegin(petscisevents,0,0,0,0);CHKERRQ(info);
    info = ISCreateComplement(this->isp, this->VSpace, &this->ispComplement);CHKERRQ(info);
    info=PetscLogEventEnd(petscisevents,0,0,0,0);CHKERRQ(info);
  }
  *isppp = this->ispComplement;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecWhichEqual"
int VecWhichEqual(Vec Vec1, Vec Vec2, IS * S)
{
  /* 
     Create an index set containing the indices of
     the vectors Vec1 and Vec2 with identical elements.
  */
  int    info;
  PetscInt i,n_same=0;
  PetscInt n,low,high,low2,high2;
  PetscInt    *same;
  PetscScalar *v1,*v2;
  MPI_Comm comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Vec1,VEC_COOKIE,1); 
  PetscValidHeaderSpecific(Vec2,VEC_COOKIE,2); 
  PetscCheckSameComm(Vec1,1,Vec2,2);


  info = VecGetOwnershipRange(Vec1, &low, &high); CHKERRQ(info);
  info = VecGetOwnershipRange(Vec2, &low2, &high2); CHKERRQ(info);

  if ( low != low2 || high != high2 )
    SETERRQ(1,"Vectors must be identically loaded over processors");

  info = VecGetLocalSize(Vec1,&n); CHKERRQ(info);

  if (n>0){
    
    if (Vec1 == Vec2){
      info = VecGetArray(Vec1,&v1); CHKERRQ(info);
      v2=v1;
    } else {
      info = VecGetArray(Vec1,&v1); CHKERRQ(info);
      info = VecGetArray(Vec2,&v2); CHKERRQ(info);
    }

    info = PetscMalloc( n*sizeof(PetscInt),&same ); CHKERRQ(info);
    
    for (i=0; i<n; i++){
      if (v1[i] == v2[i]) {same[n_same]=low+i; n_same++;}
    }
    
    if (Vec1 == Vec2){
      info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
    } else {
      info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
      info = VecRestoreArray(Vec2,&v2); CHKERRQ(info);
    }

  } else {

    n_same = 0; same=NULL;

  }

  info = PetscObjectGetComm((PetscObject)Vec1,&comm);CHKERRQ(info);
  info = ISCreateGeneral(comm,n_same,same,S);CHKERRQ(info);

  if (same) {
    info = PetscFree(same); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecWhichLessThan"
int VecWhichLessThan(Vec Vec1, Vec Vec2, IS * S)
{
  int info;
  PetscInt i;
  PetscInt n,low,high,low2,high2,n_lt=0;
  PetscInt *lt;
  PetscScalar *v1,*v2;
  MPI_Comm comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Vec1,VEC_COOKIE,1); 
  PetscValidHeaderSpecific(Vec2,VEC_COOKIE,2); 
  PetscCheckSameComm(Vec1,1,Vec2,2);

  info = VecGetOwnershipRange(Vec1, &low, &high); CHKERRQ(info);
  info = VecGetOwnershipRange(Vec2, &low2, &high2); CHKERRQ(info);

  if ( low != low2 || high != high2 )
    SETERRQ(1,"Vectors must be identically loaded over processors");

  info = VecGetLocalSize(Vec1,&n); CHKERRQ(info);

  if (n>0){

    if (Vec1 == Vec2){
      info = VecGetArray(Vec1,&v1); CHKERRQ(info);
      v2=v1;
    } else {
      info = VecGetArray(Vec1,&v1); CHKERRQ(info);
      info = VecGetArray(Vec2,&v2); CHKERRQ(info);
    }
    info = PetscMalloc( n*sizeof(PetscInt),&lt ); CHKERRQ(info);
    
    for (i=0; i<n; i++){
      if (v1[i] < v2[i]) {lt[n_lt]=high+i; n_lt++;}
    }

    if (Vec1 == Vec2){
      info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
    } else {
      info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
      info = VecRestoreArray(Vec2,&v2); CHKERRQ(info);
    }
      
  } else {
    n_lt=0; lt=NULL;
  }

  info = PetscObjectGetComm((PetscObject)Vec1,&comm);CHKERRQ(info);
  info = ISCreateGeneral(comm,n_lt,lt,S);CHKERRQ(info);

  if (lt) {
    info = PetscFree(lt); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecWhichGreaterThan"
int VecWhichGreaterThan(Vec Vec1, Vec Vec2, IS * S)
{
  int    info;
  PetscInt n,low,high,low2,high2,n_gt=0,i;
  PetscInt    *gt=NULL;
  PetscScalar *v1,*v2;
  MPI_Comm comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Vec1,VEC_COOKIE,1); 
  PetscValidHeaderSpecific(Vec2,VEC_COOKIE,2); 
  PetscCheckSameComm(Vec1,1,Vec2,2);

  info = VecGetOwnershipRange(Vec1, &low, &high); CHKERRQ(info);
  info = VecGetOwnershipRange(Vec2, &low2, &high2); CHKERRQ(info);

  if ( low != low2 || high != high2 )
    SETERRQ(1,"Vectors must be identically loaded over processors");

  info = VecGetLocalSize(Vec1,&n); CHKERRQ(info);

  if (n>0){

    if (Vec1 == Vec2){
      info = VecGetArray(Vec1,&v1); CHKERRQ(info);
      v2=v1;
    } else {
      info = VecGetArray(Vec1,&v1); CHKERRQ(info);
      info = VecGetArray(Vec2,&v2); CHKERRQ(info);
    }    

    info = PetscMalloc( n*sizeof(PetscInt), &gt ); CHKERRQ(info);
    
    for (i=0; i<n; i++){
      if (v1[i] > v2[i]) {gt[n_gt]=low+i; n_gt++;}
    }

    if (Vec1 == Vec2){
      info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
    } else {
      info = VecRestoreArray(Vec1,&v1); CHKERRQ(info);
      info = VecRestoreArray(Vec2,&v2); CHKERRQ(info);
    }
    
  } else{
    
    n_gt=0; gt=NULL;

  }

  info = PetscObjectGetComm((PetscObject)Vec1,&comm);CHKERRQ(info);
  info = ISCreateGeneral(comm,n_gt,gt,S);CHKERRQ(info);

  if (gt) {
    info = PetscFree(gt); CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecWhichBetween"
int VecWhichBetween(Vec VecLow, Vec V, Vec VecHigh, IS * S)
{
  /* 
     Creates an index set with the indices of V whose 
     elements are stictly between the corresponding elements 
     of the vector VecLow and the Vector VecHigh
  */
  int info;
  PetscInt n,low,high,low2,high2,low3,high3,n_vm=0;
  PetscInt *vm,i;
  PetscScalar *v1,*v2,*vmiddle;
  MPI_Comm comm;

  PetscValidHeaderSpecific(V,VEC_COOKIE,2); 
  PetscCheckSameComm(V,2,VecLow,1); PetscCheckSameComm(V,2,VecHigh,3);

  info = VecGetOwnershipRange(VecLow, &low, &high); CHKERRQ(info);
  info = VecGetOwnershipRange(VecHigh, &low2, &high2); CHKERRQ(info);
  info = VecGetOwnershipRange(V, &low3, &high3); CHKERRQ(info);

  if ( low!=low2 || high!=high2 || low!=low3 || high!=high3 )
    SETERRQ(1,"Vectors must be identically loaded over processors");

  info = VecGetLocalSize(VecLow,&n); CHKERRQ(info);

  if (n>0){

    info = VecGetArray(VecLow,&v1); CHKERRQ(info);
    if (VecLow != VecHigh){
      info = VecGetArray(VecHigh,&v2); CHKERRQ(info);
    } else {
      v2=v1;
    }
    if ( V != VecLow && V != VecHigh){
      info = VecGetArray(V,&vmiddle); CHKERRQ(info);
    } else if ( V==VecLow ){
      vmiddle=v1;
    } else {
      vmiddle =v2;
    }

    info = PetscMalloc( n*sizeof(PetscInt), &vm ); CHKERRQ(info);
    
    for (i=0; i<n; i++){
      if (v1[i] < vmiddle[i] && vmiddle[i] < v2[i]) {vm[n_vm]=low+i; n_vm++;}
    }

    info = VecRestoreArray(VecLow,&v1); CHKERRQ(info);
    if (VecLow != VecHigh){
      info = VecRestoreArray(VecHigh,&v2); CHKERRQ(info);
    }
    if ( V != VecLow && V != VecHigh){
      info = VecRestoreArray(V,&vmiddle); CHKERRQ(info);
    }

  } else {

    n_vm=0; vm=NULL;

  }

  info = PetscObjectGetComm((PetscObject)V,&comm);CHKERRQ(info);
  info = ISCreateGeneral(comm,n_vm,vm,S);CHKERRQ(info);

  if (vm) {
    info = PetscFree(vm); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VecWhichBetweenOrEqual"
int VecWhichBetweenOrEqual(Vec VecLow, Vec V, Vec VecHigh, IS * S)
{
  /* 
     Creates an index set with the indices of V whose 
     elements are stictly between the corresponding elements 
     of the vector VecLow and the Vector VecHigh
  */
  int info;
  PetscInt n,low,high,low2,high2,low3,high3,n_vm=0,i;
  PetscInt *vm;
  PetscScalar *v1,*v2,*vmiddle;
  MPI_Comm comm;

  PetscValidHeaderSpecific(V,VEC_COOKIE,2); 
  PetscCheckSameComm(V,2,VecLow,1); PetscCheckSameComm(V,2,VecHigh,3);

  info = VecGetOwnershipRange(VecLow, &low, &high); CHKERRQ(info);
  info = VecGetOwnershipRange(VecHigh, &low2, &high2); CHKERRQ(info);
  info = VecGetOwnershipRange(V, &low3, &high3); CHKERRQ(info);

  if ( low!=low2 || high!=high2 || low!=low3 || high!=high3 )
    SETERRQ(1,"Vectors must be identically loaded over processors");

  info = VecGetLocalSize(VecLow,&n); CHKERRQ(info);

  if (n>0){

    info = VecGetArray(VecLow,&v1); CHKERRQ(info);
    if (VecLow != VecHigh){
      info = VecGetArray(VecHigh,&v2); CHKERRQ(info);
    } else {
      v2=v1;
    }
    if ( V != VecLow && V != VecHigh){
      info = VecGetArray(V,&vmiddle); CHKERRQ(info);
    } else if ( V==VecLow ){
      vmiddle=v1;
    } else {
      vmiddle =v2;
    }

    info = PetscMalloc( n*sizeof(PetscInt), &vm ); CHKERRQ(info);
    
    for (i=0; i<n; i++){
      if (v1[i] <= vmiddle[i] && vmiddle[i] <= v2[i]) {vm[n_vm]=low+i; n_vm++;}
    }

    info = VecRestoreArray(VecLow,&v1); CHKERRQ(info);
    if (VecLow != VecHigh){
      info = VecRestoreArray(VecHigh,&v2); CHKERRQ(info);
    }
    if ( V != VecLow && V != VecHigh){
      info = VecRestoreArray(V,&vmiddle); CHKERRQ(info);
    }

  } else {

    n_vm=0; vm=NULL;

  }

  info = PetscObjectGetComm((PetscObject)V,&comm);CHKERRQ(info);
  info = ISCreateGeneral(comm,n_vm,vm,S);CHKERRQ(info);

  if (vm) {
    info = PetscFree(vm); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "ISCreateComplement"
/*@C
   ISCreateComplement - Creates the complement of the the index set

   Input Parameter:
+  S -  a PETSc IS
-  V - the reference vector space

   Output Parameter:
.  T -  the complement of S


.seealso ISCreateGeneral()

   Level: advanced
@*/
int ISCreateComplement(IS S, Vec V, IS *T){
  int info;
  PetscInt i,nis,nloc,high,low,n=0;
  const PetscInt *s;
  PetscInt *tt,*ss;
  MPI_Comm comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(S,IS_COOKIE,1); 
  PetscValidHeaderSpecific(V,VEC_COOKIE,2); 

  info = VecGetOwnershipRange(V,&low,&high); CHKERRQ(info);
  info = VecGetLocalSize(V,&nloc); CHKERRQ(info);
  info = ISGetLocalSize(S,&nis); CHKERRQ(info);
  info = ISGetIndices(S, &s); CHKERRQ(info);
  info = PetscMalloc( nloc*sizeof(PetscInt),&tt ); CHKERRQ(info);
  info = PetscMalloc( nloc*sizeof(PetscInt),&ss ); CHKERRQ(info);

  for (i=low; i<high; i++){ tt[i-low]=i; }

  for (i=0; i<nis; i++){ tt[s[i]-low] = -2; }
  
  for (i=0; i<nloc; i++){
    if (tt[i]>-1){ ss[n]=tt[i]; n++; }
  }

  info = ISRestoreIndices(S, &s); CHKERRQ(info);
  
  info = PetscObjectGetComm((PetscObject)S,&comm);CHKERRQ(info);
  info = ISCreateGeneral(comm,n,ss,T);CHKERRQ(info);
  
  if (tt) {
    info = PetscFree(tt); CHKERRQ(info);
  }
  if (ss) {
    info = PetscFree(ss); CHKERRQ(info);
  }

  PetscFunctionReturn(0);
}
#endif

