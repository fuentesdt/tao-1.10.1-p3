#include "esi/ESI.h"
// #include "tao_general.h"
// #include "/home/benson/Trilinos/src/petra_esi/Petra_ESI.h"

#define TaoFunctionBegin info=1
#define TaoFunctionReturn(a) return a;
#define CHKERRQ(n) \
  {\
  if (n) {return n;} \
  }

#define SETERRQ(a,b) \
  {\
  {return 1;} \
  }
 

#undef __FUNCT__
#define __FUNCT__ "ESIshiftDiagonal"
int ESIshiftDiagonal( esi::MatrixRowWriteAccess<double,int>* esimatwrite, double c){
  esi::ErrorCode info;
  double cc=c;
  esi::IndexSpace<int>* rowmap = NULL, *colmap = NULL;
  int i,row,localSize,firstLocalEqn;

  TaoFunctionBegin;
  info = esimatwrite->getIndexSpaces(rowmap, colmap); CHKERRQ(info);

  info = rowmap->getLocalSize(localSize); CHKERRQ(info);
  info = rowmap->getLocalPartitionOffset(firstLocalEqn); CHKERRQ(info);
  for (i=0;i<localSize; i++){
    row = firstLocalEqn+i;
    cc=c;
    info = esimatwrite->sumIntoRow(row, &cc, &row, 1); CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ESIaddDiagonal"
int ESIaddDiagonal(esi::Object *esiobj2,  esi::Vector<double,int> *X){
  esi::ErrorCode info;
  double cc;
  esi::IndexSpace<int>* rowmap = NULL, *colmap = NULL;
  esi::MatrixRowWriteAccess<double,int>* esimatwrite;
  int i,row,localSize,firstLocalEqn;
  double *x;

  TaoFunctionBegin;
  esimatwrite = dynamic_cast<esi::MatrixRowWriteAccess<double,int>*>(esiobj2);
  if (esimatwrite==NULL){
    SETERRQ(56,"Operation not defined");
  } else {
    info = esimatwrite->getIndexSpaces(rowmap, colmap); CHKERRQ(info);
    rowmap->getLocalSize(localSize);
    rowmap->getLocalPartitionOffset(firstLocalEqn);

    info=X->getCoefPtrReadWriteLock(x); CHKERRQ(info);
    for (i=0;i<localSize; i++){
      row = firstLocalEqn+i; cc=x[i];
      info = esimatwrite->sumIntoRow(row, &cc, &row, 1); CHKERRQ(info);
    }
      info=X->releaseCoefPtrLock(x); CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ESIsetDiagonal"
int ESIsetDiagonal(esi::Object *esiobj2,  esi::Vector<double,int> *X){
  esi::ErrorCode info;
  double cc;
  esi::MatrixRowWriteAccess<double,int>* esimatwrite;
  esi::IndexSpace<int>* rowmap = NULL, *colmap = NULL;
  int i,row,localSize,firstLocalEqn;
  double *x;

  TaoFunctionBegin;
  esimatwrite = dynamic_cast<esi::MatrixRowWriteAccess<double,int>*>(esiobj2);
  if (esimatwrite==NULL){
    SETERRQ(56,"Operation not defined");
  } else {
    info = esimatwrite->getIndexSpaces(rowmap, colmap); CHKERRQ(info);
    info = rowmap->getLocalSize(localSize); CHKERRQ(info);
    info = rowmap->getLocalPartitionOffset(firstLocalEqn); CHKERRQ(info);

    info=X->getCoefPtrReadWriteLock(x); CHKERRQ(info);
    for (i=0;i<localSize; i++){
      row = firstLocalEqn+i; cc=x[i];
      info = esimatwrite->copyIntoRow(row, &cc, &row, 1); CHKERRQ(info);
    }
    info=X->releaseCoefPtrLock(x); CHKERRQ(info);
  }

  TaoFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ESIgetDiagonal"
int ESIgetDiagonal(esi::Object *esiobj2,  esi::Vector<double,int> *X){
  esi::ErrorCode info;
  double dtmp;
  esi::MatrixRowReadAccess<double,int>* esimatread;
  esi::IndexSpace<int>* rowmap = NULL, *colmap = NULL;
  int i,j,m,n,row,rowLen,localSize,firstLocalEqn;
  double *x, *coefs;
  int *colIndices;

  TaoFunctionBegin;
  esimatread = dynamic_cast<esi::MatrixRowReadAccess<double,int>*>(esiobj2);
  if (esimatread==NULL){
    SETERRQ(56,"Operation not defined");
  } else {
    
    info = esimatread->getIndexSpaces(rowmap, colmap); CHKERRQ(info);
    info = rowmap->getLocalSize(localSize); CHKERRQ(info);
    info = rowmap->getLocalPartitionOffset(firstLocalEqn); CHKERRQ(info);
    
    info=esimatread->getGlobalSizes(m,n); CHKERRQ(info);
    coefs=new double[m+1];
    colIndices = new int[m+1];
    info=X->getCoefPtrReadWriteLock(x); CHKERRQ(info);
    for (i=0;i<localSize; i++){
      row = firstLocalEqn+i;
      dtmp=0.0;
      info=esimatread->getRowNonzeros(row, rowLen); CHKERRQ(info);
      info=esimatread->copyOutRow(row,coefs,colIndices,m,rowLen); CHKERRQ(info);
      for (j=0;j<rowLen;j++){
	if (colIndices[j]==row){
	  dtmp=coefs[j]; 
	  break;
	}
      }
      x[i]=dtmp;
    }
    delete coefs;
    delete colIndices;
    info=X->releaseCoefPtrLock(x); CHKERRQ(info);    
  }

  TaoFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ESIRowScale"
int ESIRowScale( esi::Object *esiobj, esi::Vector<double,int>*X){
  esi::ErrorCode info;
  double dtmp;
  esi::IndexSpace<int>* rowmap = NULL;
  esi::MatrixRowReadAccess<double,int>* esimatread;
  esi::MatrixRowWriteAccess<double,int>* esimatwrite;
  int i,j,m,n,row,localSize,firstLocalEqn,rowLen;
  int *colIndices;
  double *x,*coefs;

  TaoFunctionBegin;
  esimatwrite = dynamic_cast<esi::MatrixRowWriteAccess<double,int>*>(esiobj);
  esimatread = dynamic_cast<esi::MatrixRowReadAccess<double,int>*>(esiobj);
  if (esimatwrite==NULL || esimatread==NULL){
    SETERRQ(56,"Operation not defined");
  } else {
    info = X->getIndexSpace(rowmap); CHKERRQ(info);
    info = rowmap->getLocalSize(localSize); CHKERRQ(info);
    info = rowmap->getLocalPartitionOffset(firstLocalEqn); CHKERRQ(info);
    info=esimatread->getGlobalSizes(m,n); CHKERRQ(info);
    coefs=new double[m+1];
    colIndices = new int[m+1];
    info=X->getCoefPtrReadWriteLock(x); CHKERRQ(info);
    for (i=0;i<localSize; i++){
      row = firstLocalEqn+i;
      dtmp=x[i];
      info=esimatread->getRowNonzeros(row, rowLen); CHKERRQ(info);
      info=esimatread->copyOutRow(row,coefs,colIndices,m,rowLen); CHKERRQ(info);
      for (j=0;j<rowLen;j++){
	  coefs[j]*=dtmp; 
      }
      info = esimatwrite->copyIntoRow(row,coefs,colIndices,rowLen); CHKERRQ(info);
    }

    delete coefs;
    delete colIndices;
    info=X->releaseCoefPtrLock(x); CHKERRQ(info);    

  }
  TaoFunctionReturn(0);
}


#include <iostream.h>
#undef __FUNCT__
#define __FUNCT__ "TaoPrint_ESI_RowMatrix"
int TaoPrint_ESI_RowMatrix(esi::Object * esiobj)
{
   esi::IndexSpace<int>* rowmap = NULL, *colmap = NULL;
   esi::MatrixRowReadAccess<double,int>* esimat=NULL;
   esi::ErrorCode info;

   TaoFunctionBegin;
   esimat= dynamic_cast<esi::MatrixRowReadAccess<double,int>*>(esiobj);
   if (esimat==NULL){
     SETERRQ(56,"Operation not defined");
   }

   info = esimat->getIndexSpaces(rowmap, colmap); CHKERRQ(info);

   //we need to know how many local equations there are, and what the first
   //local equation is. (Since we happen to know that the maps were built with
   //an indexBase of 0, then firstLocalEqn == 'getLocalPartitionOffset'.)
   int localSize, firstLocalEqn, localRank;
   info = rowmap->getLocalSize(localSize); CHKERRQ(info);
   info = rowmap->getLocalPartitionOffset(firstLocalEqn); CHKERRQ(info);
   info = rowmap->getLocalPartitionRank(localRank); CHKERRQ(info);

   int numRows, numCols;
   info = esimat->getGlobalSizes(numRows, numCols); CHKERRQ(info);
   cout << localRank << ": global rows " << numRows << ", global cols "
	<< numCols << endl;

   int* colIndices = new int[numRows+1];
   double * coefs = new double[numRows+1]; 

   for(int i=0; i<localSize; i++) {
      //first, make sure our colIndices and coefs arrays are the right length.
      //(This operation is a no-op if the array is already the right size or
      //bigger.)

      int rowLen;
      int row = firstLocalEqn+i;
      info = esimat->getRowNonzeros(row, rowLen); CHKERRQ(info);

      info = esimat->copyOutRow(row, coefs, colIndices,
				numCols, rowLen); CHKERRQ(info);
      cout << localRank << ": row " << row << ": ";
      for(int j=0; j<rowLen; j++) {
	cout << "("<<colIndices[j]<<","<<coefs[j]<<") ";
      }
      cout << endl;

   }
   delete colIndices;
   delete coefs;
   TaoFunctionReturn(0);
}
