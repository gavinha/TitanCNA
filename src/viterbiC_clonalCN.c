/** viterbiC_clonalCN.c
* C implementation of viterbi algorithm
* Special features include:
* 1) joint genotype and clonal states
* 2) specific genotype and clonal transition probabilities
*
* This function WORKS WITH LOGARITHMS.
* 
* author: Gavin Ha <gha@bccrc.ca>
*          Dept of Molecular Oncolgy
*          British Columbia Cancer Agency
*          University of British Columbia
* date  : November 9, 2012
*
* 
**/

#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>


void copyVectorC(double *, double *, unsigned int);
void addVectors(double *, double *, double *, unsigned int);
void setVectorToValue(double *, double, unsigned int);
void setVectorToValue_int(int *, int, unsigned int);
void maxVectorInPlace(double *, int *, double *, unsigned int);
void initializeTxnV(double * transSlice, unsigned int K);
void logMatrixInPlace(double * A, unsigned int K);
void outputMatrixV(double * A, unsigned int K);
double distanceTransitionFunctionV(double, double, double);
void preparePositionSpecificMatrix(double *, unsigned int, unsigned int, double *, double *, double, double, unsigned int, unsigned int);

SEXP viterbiC_clonalCN(SEXP piGiZi, SEXP py, SEXP copyNumKey, SEXP zygosityKey, SEXP numClust, SEXP positions, SEXP zStrength, SEXP txnLen, SEXP useOutlier) {
  double * prior, * obslik, * transSlice;
  int K, Z, T, numUnitStates;
  double * CT, * ZS, * posn, * txnExpLen, *txnZstrength, *numClusters, * outlier;
  double * delta;
  int * psi;
  int t, j;
  double rhoG, rhoZ;
  double * d; /* buffer */
  
  PROTECT(piGiZi = AS_NUMERIC(piGiZi));  // K*Z-by-1 
  PROTECT(py = AS_NUMERIC(py));    // K*Z-by-T, 2D-array obslik
  PROTECT(copyNumKey = AS_NUMERIC(copyNumKey)); // K-by-1
  PROTECT(zygosityKey = AS_NUMERIC(zygosityKey)); // K-by-1
  PROTECT(numClust = AS_NUMERIC(numClust)); // scalar
  PROTECT(positions = AS_NUMERIC(positions)); // T-by-1
  PROTECT(zStrength = AS_NUMERIC(zStrength)); // scalar
  PROTECT(txnLen = AS_NUMERIC(txnLen)); // scalar
	PROTECT(useOutlier = AS_NUMERIC(useOutlier)); //scalar
  prior = NUMERIC_POINTER(piGiZi); 
  obslik = NUMERIC_POINTER(py); 
  CT = NUMERIC_POINTER(copyNumKey);
  ZS = NUMERIC_POINTER(zygosityKey);
  numClusters = NUMERIC_POINTER(numClust);
  posn = NUMERIC_POINTER(positions); 
  txnZstrength = NUMERIC_POINTER(zStrength); 
  txnExpLen = NUMERIC_POINTER(txnLen); 
	outlier = NUMERIC_POINTER(useOutlier); 
    
  /* Check size of initial state distribution */
  K = GET_LENGTH(piGiZi);
  T = GET_LENGTH(positions);
  Z = (int)numClusters[0];
  numUnitStates = (int)K/Z;    
    
  if (INTEGER(GET_DIM(py))[0] != K || INTEGER(GET_DIM(py))[1] != T){
    error("viterbiC_clonalCN: The obslik must be %d-by-%d dimension.",K,T);
  }
/*
  if (GET_LENGTH(copyNumKey) != K/Z){
    error("viterbiC_clonalCN: The copy number vector must be length %d.",numUnitStates);
  }
*/  
  if (GET_LENGTH(positions) != T){
    error("viterbiC_clonaCN: The positions vector must be of size %d-by-1.",T);
  }
  
  delta = malloc(K*T*sizeof(double));
  psi = malloc(K*T*sizeof(double));
  // setup the only variable we wish to return, optimal path
  SEXP path_data;
  PROTECT(path_data = allocVector(INTSXP,T));
  int * path;
  path = INTEGER_POINTER(path_data);
  
  t = 0;
  addVectors(delta + t*K, prior, obslik + t*K, K);
  setVectorToValue_int(psi + t*K, 0, K);
  
  d = malloc(K*sizeof(double));
  transSlice = (double *)malloc(K*K*sizeof(double));
    
  /* forward */
  for(t=1;t<T;++t) { /* position */
      /* Each iteration, we overwrite transSlice with transmat
       * to start over when at a new probe */
      initializeTxnV(transSlice, K);
      /* modify transSlice inplace by adding position-specific probs */
      rhoG = 1.0 - distanceTransitionFunctionV(posn[t-1],posn[t],txnExpLen[0]);    
      rhoZ = 1.0 - distanceTransitionFunctionV(posn[t-1],posn[t],txnZstrength[0]);
      preparePositionSpecificMatrix(transSlice, K, numUnitStates, CT, ZS, rhoG, rhoZ, outlier[0], 0);   
      logMatrixInPlace(transSlice, K);        
      for(j=0;j<K;++j) { /* column */          
          addVectors(d, delta + (t-1)*K, transSlice + j*K , K);
          maxVectorInPlace(delta + j + t*K, psi + j + t*K, d, K);
          delta[j+t*K] += obslik[j+t*K];
      }
  
  }
  
  
  /* backward */
  t = T-1;
  maxVectorInPlace(d, path + t, delta + t*K, K); /* using the first value of d to store junk */
  
  for(t=T-2;t>=0;--t) {
      path[t] = psi[path[t+1] + (t+1)*K];
  }
     
  /* Be careful to add +1 to path values. This is because C starts from 0
   * and R starts from 1.
   */
  for(t=0;t<T;++t)
      path[t] = (double)(path[t]+1);
  
  free(delta); free(psi); free(d);
  UNPROTECT(10);

  return path_data;
}



/* Position specific distance used in transition matrix
 * returns double
 */
double distanceTransitionFunctionV(double prevPosn, double curPosn, double L) {
    double distance = 0;
    double rho = 0;
    distance = curPosn - prevPosn + 1.0; /* won't encounter next chr */
    rho = (1.0/2.0)*(1.0-exp(-distance/(2.0*L)));
    return rho;
}

/* Method to make a deep copy of a K-by-K matrix */
void initializeTxnV(double * transSlice, unsigned int K) {
    unsigned int i, j;
    
    /* deep copy transmat to transSlice */
    for (i=0;i<K;i++) /* rows */
        for (j=0;j<K;j++) /* columns */
            transSlice[i + j*K] = 0;
}


/* Output matrix for debugging purposes */
void outputMatrixV(double * A, unsigned int K) {
    unsigned int i, j;
    
    for (i=0;i<K;i++) { /*rows*/
        for (j=0;j<K;j++) { /*cols*/
            //printf("%f\t", A[K*j+i]);
        }
        //printf("\n");
    }
}

void copyVectorC(double * Out, double * In, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        Out[i] = In[i];
    
    return;
}

void addVectors(double * Out, double * u, double * v, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        Out[i] = u[i] + v[i];
    
    return;
    
}

void setVectorToValue(double * A, double value, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        A[i] = value;
    
    return;
}

void setVectorToValue_int(int * A, int value, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        A[i] = value;
    
    return;
}


void maxVectorInPlace(double * Out_value, int * Out_index, double * A, unsigned int L) {
    unsigned int i;
    double maxvalue;
    int index;
    
    maxvalue = A[0];
    index = 0;
    
    for(i=1;i<L;++i) {
        if (maxvalue < A[i]) {
            index = i;
            maxvalue = A[i];
        }
    }
    
    *Out_value = maxvalue;
    *Out_index = index;
    
    return;
}

