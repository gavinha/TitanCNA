/** getPositionOverlap.c
 * C implementation of finding overlap given
 * 1) Set of positions of interest
 * 2) Larger list of start and stop boundaries we wish to find overlap in
 *
 * author: Gavin Ha <gha@bccrc.ca>
 *          Dept of Molecular Oncolgy
 *          British Columbia Cancer Agency
 *          University of British Columbia
 * date  : November 7, 2012
 *
 **/

#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

SEXP getPositionOverlapC(SEXP posn, SEXP start, SEXP stop){
  double *posns, *cnStart, *cnStop; /* inputs */ 
  int T, Tcn, i, j;  
  
  PROTECT(posn = AS_NUMERIC(posn));   // vector
  PROTECT(start = AS_NUMERIC(start)); // vector
  PROTECT(stop = AS_NUMERIC(stop));  // vector
  posns = NUMERIC_POINTER(posn);
  cnStart = NUMERIC_POINTER(start);
  cnStop = NUMERIC_POINTER(stop);
  
  T = GET_LENGTH(posn);
  Tcn = GET_LENGTH(start);  
  
  SEXP outIndicesData;
  double *outIndices;
  PROTECT(outIndicesData = allocVector(REALSXP,T));
  outIndices = NUMERIC_POINTER(outIndicesData);
  
  for (i=0; i<T; i++){  // for each data point  
    outIndices[i] = 0;
    for (j=0; j<Tcn; j++){  // for each copy number segment 
      if (((int)cnStart[j]<=(int)posns[i]) && ((int)posns[i]<=(int)cnStop[j])){
          outIndices[i] = j+1;
          break;           
      }
    }
  }    
  
  UNPROTECT(4);
  return(outIndicesData);

  
}
