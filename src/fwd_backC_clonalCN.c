/** fwd_backC_clonalCN.c
 * C implementation of forwards-backwards algorithm for R-2.15
 * Special features include:
 * 1) joint genotype and clonal states
 * 2) specific genotype and clonal transition probabilities treated as mega-variable
 *    i.e. state space is K*Z where K is number of genotypes and Z is number of clonal states
 * Precondition:
 * py and piGiZi must be in log space
 *
 * author: Gavin Ha <gavinha@gmail.com>
 *          Dept of Molecular Oncolgy
 *          British Columbia Cancer Agency
 *          University of British Columbia
 * date  : June 4, 2014
 *
 **/
 
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

double normalizeInPlace(double *, unsigned int);
double logSumNormalizeInPlace(double *, unsigned int);
void multiplyInPlace(double *, double *, double *, unsigned int);
void addInPlace(double *, double *, double *, unsigned int);
void multiplyMatrixInPlace(double *, double *, double *, unsigned int);
void logSumInPlace(double *, double *, double *, unsigned int);
double logsumexp(double *, unsigned int);
void logMatrixInPlace(double *, unsigned int);
void transposeSquareInPlace(double *, double *, unsigned int);
void outerProductUVInPlace(double *, double *, double *, unsigned int);
void componentVectorMultiplyInPlace(double *, double *, double *, unsigned int);
void preparePositionSpecificMatrix(double *, unsigned int, unsigned int, double *, double *, double, double, unsigned int, unsigned int);
void initializeTxn(double *, unsigned int);
void outputMatrix(double *, unsigned int);
void outputVector(double *, unsigned int);
double distanceTransitionFunction(double, double, double);

SEXP fwd_backC_clonalCN(SEXP piGiZi, SEXP py, SEXP copyNumKey, SEXP zygosityKey, SEXP numClust, SEXP positions, SEXP zStrength, SEXP txnLen, SEXP useOutlier) {
  double * init_state_distrib, * obslik;
  double * CT, * ZS, * posn, * txnExpLen, * txnZstrength, * numClusters, * outlier;
  int K, T, Z, numUnitStates; 
  int t, d;
  double rhoG = 0.0;
  double rhoZ = 0.0;
    
  PROTECT(piGiZi = AS_NUMERIC(piGiZi));  // K*Z-by-1 
  PROTECT(py = AS_NUMERIC(py));    // K*Z-by-T, 2D-array obslik
  PROTECT(copyNumKey = AS_NUMERIC(copyNumKey)); // K-by-1
  PROTECT(zygosityKey = AS_NUMERIC(zygosityKey)); // K-by-1
  PROTECT(numClust = AS_NUMERIC(numClust)); // scalar
  PROTECT(positions = AS_NUMERIC(positions)); // T-by-1
  PROTECT(zStrength = AS_NUMERIC(zStrength)); // scalar
  PROTECT(txnLen = AS_NUMERIC(txnLen)); // scalar
	PROTECT(useOutlier = AS_NUMERIC(useOutlier)); //scalar
  init_state_distrib = NUMERIC_POINTER(piGiZi);
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
	numUnitStates = (int)(K)/Z;
 
  if (INTEGER(GET_DIM(py))[0] != K || INTEGER(GET_DIM(py))[1] != T){
    error("fwd_backC_clonalCN: The obslik must be %d-by-%d dimension.",K,T);
  }
/*
  if (GET_LENGTH(copyNumKey) != numUnitStates){
    error("fwd_backC_clonalCN: The copy number vector must be length %d.",numUnitStates);
  } 
*/  
  if (GET_LENGTH(positions) != T){
    error("fwd_backC_clonaCN: The positions vector must be of size %d-by-1.",T);
  }
  
  //SEXP scale_data, alpha_data, beta_data;
  SEXP gamma_data, loglik_data;
  double * scale, * alpha, * beta, * gamma, * loglik, * m;
  //PROTECT(scale_data = NEW_NUMERIC(T));
  //PROTECT(alpha_data = allocMatrix(REALSXP, K, T));
  //PROTECT(beta_data = allocMatrix(REALSXP, K, T));
  PROTECT(gamma_data = allocMatrix(REALSXP, K, T));
  PROTECT(loglik_data = NEW_NUMERIC(1));
  //scale = NUMERIC_POINTER(scale_data);
  //alpha = NUMERIC_POINTER(alpha_data);
  //beta = NUMERIC_POINTER(beta_data);
  gamma = NUMERIC_POINTER(gamma_data);
  loglik = NUMERIC_POINTER(loglik_data);
    
  scale = malloc(T*sizeof(double));    
  alpha = malloc(K*T*sizeof(double));
  beta = malloc(K*T*sizeof(double));
  //gamma = malloc(K*T*sizeof(double));
    
  /* Use transSlice as the temporary txn matrix to modify
   * with position-specific; each iteration, we overwrite it with transmat
   * to start over at a new probe */
  double * transmatT, * transSlice;  
  transSlice = malloc(K*K*sizeof(double));    
  transmatT = malloc(K*K*sizeof(double));    
    
    /********* Forward. ********/
   t = 0;
   /*multiplyInPlace(alpha + t*K, init_state_distrib, obslik + t*K, K);*/
   addInPlace(alpha + t*K, init_state_distrib, obslik + t*K, K);
   scale[t] = logSumNormalizeInPlace(alpha + t*K, K);
    
   m = malloc(K*sizeof(double));

  for(t=1;t<T;++t){
  /* Each iteration, we overwrite transSlice with transmat
  * to start over when at a new probe */
    initializeTxn(transSlice, K);     
    /* modify transSlice inplace by adding position-specific probs */
    rhoG = 1.0 - distanceTransitionFunction(posn[t-1],posn[t],txnExpLen[0]);
    rhoZ = 1.0 - distanceTransitionFunction(posn[t-1],posn[t],txnZstrength[0]);        
    preparePositionSpecificMatrix(transSlice, K, numUnitStates, CT, ZS, rhoG, rhoZ, outlier[0], 0);   
    transposeSquareInPlace(transmatT, transSlice, K);
    logMatrixInPlace(transmatT, K);
    /*multiplyMatrixInPlace(m, transmatT, alpha + (t-1)*K, K);*/
    logSumInPlace(m, transmatT, alpha + (t-1)*K, K);
    /*multiplyInPlace(alpha + t*K, m, obslik + t*K, K); */
    addInPlace(alpha + t*K, m, obslik + t*K, K);
    //printf("Forward t=%d\n",t);
    //outputVector(alpha + t*K, K);
    scale[t] = logSumNormalizeInPlace(alpha + t*K, K); 
		//printf("Forward t=%d\tScale=%0.2f\n",t,scale[t]);
    //outputVector(alpha + t*K, K);
  }
  
  loglik[0] = 0;
  for(t=0;t<T;++t)
  	loglik[0] += scale[t]; /* Already in log space */
  	/*loglik[0] += log(scale[t]);*/

  /********* Backward. ********/
  
  
  t = T-1;
  /* I don't think we need to initialize beta to all zeros. */
  for(d=0;d<K;++d) {
    beta[d + t*K] = 0;
    gamma[d + t*K] = alpha[d + t*K];
  }
  
  double * b;
  b = malloc(K*sizeof(double));
  
  for(t=(T-2);t>=0;--t) {        
    /* setting beta */
    /*multiplyInPlace(b, beta + (t+1)*K, obslik + (t+1)*K, K);*/
    addInPlace(b, beta + (t+1)*K, obslik + (t+1)*K, K);
    /* Using "m" again instead of defining a new temporary variable.
     * We using a lot of lines to say
     * beta(:,t) = normalize(transmat * b);
     */
    initializeTxn(transSlice, K);
    /* modify transSlice inplace by adding position-specific probs */
    rhoG = 1.0 - distanceTransitionFunction(posn[t],posn[t+1],txnExpLen[0]);
    rhoZ = 1.0 - distanceTransitionFunction(posn[t],posn[t+1],txnZstrength[0]);
    preparePositionSpecificMatrix(transSlice, K, numUnitStates, CT, ZS, rhoG, rhoZ, outlier[0], 0);        
    logMatrixInPlace(transSlice, K);
    /*multiplyMatrixInPlace(m, transSlice, b, K);*/
    logSumInPlace(m, transSlice, b, K);
    //printf("Backward t=%d\n",t);
    //outputVector(m, K);
    logSumNormalizeInPlace(m, K);
		//printf("Backward t=%d\tScale=%0.2f\n",t,sumBack);
    //outputVector(m, K);
    
    for(d=0;d<K;++d) { beta[d + t*K] = m[d]; }
    /* using "m" again as valueholder */   
    /* setting gamma */        
    /*multiplyInPlace(m, alpha + t*K, beta + t*K, K);*/
    addInPlace(m, alpha + t*K, beta + t*K, K);
    logSumNormalizeInPlace(m, K);
    //printf("FWBK t=%d\tScale=%0.2f\n",t,sumBack);
    //outputVector(m, K);
    
    for(d=0;d<K;++d) { gamma[d + t*K] = m[d]; } 
  	//printf("GAMMA t=%d\n",t);
    //outputVector(gamma + t*K, K);
  }
      
  free(b); free(m); 
  free(scale); free(transmatT); free(transSlice);
  free(alpha); free(beta);     
    
  SEXP list, list_names;
  char *names[2] = {"rho", "loglik"};
  PROTECT(list_names = allocVector(STRSXP, 2));
  for (int i = 0; i < 2; ++i) {
    SET_STRING_ELT(list_names, i, mkChar(names[i]));
  }
  PROTECT(list = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(list, 0, gamma_data); 
  SET_VECTOR_ELT(list, 1, loglik_data);
  setAttrib(list, R_NamesSymbol, list_names);
  UNPROTECT(13);
  return list;
  
}

/* Method to assign transSlice as the position specific matrix
 preparePositionSpecificMatrix(double * transSlice, unsigned int, double * ct, double * ZS, double, unsigned int) */
void preparePositionSpecificMatrix(double * transSlice, unsigned int K, unsigned int numUnitStates, double * CT, double * ZS, double rhoG, double rhoZ, unsigned int OUTLIERSTATE, unsigned int boolTest) {
    unsigned int i, j, unitI, unitJ;
    double sum, z1, z2, iZS, jZS;   
    //Z = (int)K/numUnitStates;
    /* Add the distance to our output matrix, in place 
     * Also multiple by copy number values */       
    for (i = 0; i < K; i++){ /* rows */
        //lots of code to figure out what the zygosity status of the states are based on whethe we 
        // are using the outlier state
        if (OUTLIERSTATE==1){
					if (i==0){//garbage state 
						z1 = 0;  unitI = -1;
						iZS = -100;
					}else{
						z1 = ceil(((double)i)/numUnitStates); /* cluster number */
						unitI = (int)(i-1)%(numUnitStates); /* unit state */
						iZS = ZS[unitI];
					}
        }else{
					z1 = ceil(((double)i+1)/numUnitStates); /* cluster number */
					unitI = (int)i%(numUnitStates);  /* unit state */
					iZS = ZS[unitI];
				} 
        for (j = 0; j < K; j++){ /* columns */
          	if (OUTLIERSTATE==1){
							if (j==0){//garbage state 
								z2 = 0;  unitJ = -1;
								jZS = -100;
							}else{
								z2 = ceil(((double)j)/numUnitStates); /* cluster number */
								unitJ = (int)(j-1)%(numUnitStates); /* unit state */
								jZS = ZS[unitJ];
							}
        		}else{
							z2 = ceil(((double)j+1)/numUnitStates); /* cluster number */
							unitJ = (int)j%(numUnitStates); /* unit state */
							jZS = ZS[unitJ];
						} 
            //printf("i=%d\tj=%d\tz1=%f\tz2=%f\tZS[%d]=%f\tZS[%d]=%f\n",i,j,z1,z2,unitI,iZS,unitJ,jZS);
						//transitions to same state or same zygosity status
						
            /** GENOTYPE TRANSITION **/
            if (iZS==jZS){
                transSlice[i + j*K] = rhoG; 
            }else{
	              transSlice[i + j*K] = (1.0-rhoG)/((double)K-1.0); 
            }

			/** CLONAL CLUSTER TRANSITION **/
			//if (K > numUnitStates){ //only use if more than 1 cluster
	    	//same clust or het (-1) 
				if(z1 == z2 || jZS == -1){ // || CT[unitJ] >= 4){  
					transSlice[i + j*K] = transSlice[i + j*K] * rhoZ;
				}else{  //different cluster (except for diploid HET)
					transSlice[i + j*K] = transSlice[i + j*K] * (1.0-rhoZ);
				}
			//}
        }
    }          
   
      
    if (boolTest) {
        //printf("fwd_backC: RhoG:\t%f\tRhoZ:\t%f\tOutlier: %d\tNumUnitStates=%d\n",rhoG,rhoZ,OUTLIERSTATE,numUnitStates);
        //printf("fwd_backC: Raw Matrix: \n");        
        //outputMatrix(transSlice, K);
    }
    
    /* Normalize matrix by rows */        
    for (i=0;i<K;i++) { /* rows */
        sum = 0;
        for (j=0;j<K;j++) { /* columns */
            sum += transSlice[i+j*K];
        }        
        if (sum > 0){
            for (j=0;j<K;j++) { /* columns */            
                transSlice[i+j*K] /= sum;
            }
        }
        
    }
    
    if (boolTest) {
        //printf("fwd_backC: Normalized Matrix: \n");
        //outputMatrix(transSlice, K);
    }
    
}


/*
 * Create a deep copy of a matrix
 */
void initializeTxn(double * transSlice, unsigned int K) {
    unsigned int i, j;
    
    /* deep copy transmat to transSlice */
    for (i=0;i<K;i++) /* rows */
        for (j=0;j<K;j++) /* columns */
	        transSlice[i + j*K] = 0;
}

/* Output matrix for debugging purposes */
void outputMatrix(double * A, unsigned int K) {
    unsigned int i, j;
    
    for (i=0;i<K;i++) { /*rows*/
        for (j=0;j<K;j++) { /*cols*/
            //printf("%0.2e\t", A[K*j+i]);
        }
        //printf("\n");
    }
    
}

/* Output vector for debugging purposes */
void outputVector(double * A, unsigned int K) {
    unsigned int i;
    
    for (i=0;i<K;i++) { /*rows*/        
            //printf("%0.2f\t", A[i]);
    }
    //printf("\n");
}

/* Position specific distance used in transition matrix
 * returns double
 */
double distanceTransitionFunction(double prevPosn, double curPosn, double L) {
    double distance = 0;
    double rho = 0;
    distance = curPosn - prevPosn + 1.0; /* won't encounter next chr */
    rho = (1.0/2.0)*(1.0-exp(-distance/(2.0*L)));
    return rho;
}

/* Returns the normalization constant used.
 */
double normalizeInPlace(double * A, unsigned int N) {
    unsigned int n;
    double sum = 0;
    
    for(n=0;n<N;++n) {
        sum += A[n];
        if (A[n] < 0) {
	          //printf("OurNegVal=%0.2e\n",A[n]);
            error("We don't want to normalize if A contains a negative value. This is a logical error.");
        }
    }    
    if (sum == 0)
        error("We are asked to normalize a section of a vector containing only zeros.");
    else {
        for(n=0;n<N;++n)
            A[n] /= sum;
    }
    return sum;
}

/* Returns the normalization constant used.
 * Input vector is in log space; use logsumexp function
 */
double logSumNormalizeInPlace(double * A, unsigned int N) {
	unsigned int i;
	double sum = 0;
	
	sum = logsumexp(A,N);
	for (i=0; i<N; i++){
		A[i] -= sum;
	}
	
	return sum;
}

void multiplyInPlace(double * result, double * u, double * v, unsigned int K) {
    unsigned int n;
    
    for(n=0;n<K;++n){
        result[n] = u[n] * v[n];
    }
    return;
}

void addInPlace(double * result, double * u, double * v, unsigned int K) {
    unsigned int n;
    
    for(n=0;n<K;++n){
        result[n] = u[n] + v[n];
    }
    return;
}

void multiplyMatrixInPlace(double * result, double * trans, double * v, unsigned int K) {
    
    unsigned int i, d;
    
    for(d=0;d<K;++d) {
        result[d] = 0;
        for (i=0;i<K;++i){
            result[d] += trans[d + i*K] * v[i];	    
        }
    }
    return;
}

void logSumInPlace(double * result, double * trans, double * v, unsigned int K) {
    
    unsigned int i, d;
    
    for(d=0;d<K;++d) {
        result[d] = 0;
        double * sums = malloc(K*sizeof(double)); /* keep track of trans+v */
        for (i=0;i<K;++i){
            sums[i] = trans[d + i*K] + v[i];
        }
        result[d] = logsumexp(sums, K);
        free(sums);
    }
    
    return;
}

/** 
	Computes the log of the sum of exponentials 
	log(sum(exp(x))) = A + log(sum(exp(x-A))); where A=max(x) and x is vector to be summed
	Code reference: http://stackoverflow.com/questions/4169981/logsumexp-implementation-in-c
**/
double logsumexp(double *x, unsigned int K) {
	double max_x = x[0];
	double sum = 0.0;
	unsigned int i;
	
	/* Find max value of x*/
	for (i=0; i<K; i++){
		if (x[i] > max_x){
			max_x = x[i];
		}
	}
	/* compute sum(exp(x-max_x)) */
	for (i=0; i<K; i++){
		sum += exp(x[i] - max_x);
	}	
	/* return max_x+log(sum(exp(x-max_x))) */
	return max_x + log(sum);
}

/* logs each element in a K-by-K matrix, A */
void logMatrixInPlace(double * A, unsigned int K) {
    unsigned int i, j;
    for (i=0;i<K;i++) /* rows */
        for (j=0;j<K;j++) /* columns */
            A[i + j*K] = log(A[i + j*K]);
}

void transposeSquareInPlace(double * out, double * in, unsigned int K) {
    
    unsigned int i, j;
    
    for(i=0;i<K;++i){
        for(j=0;j<K;++j){
            out[j+i*K] = in[i+j*K];
        }
    }
    return;
}

void outerProductUVInPlace(double * Out, double * u, double * v, unsigned int K) {
    unsigned int i, j;
    
    for(i=0;i<K;++i){
        for(j=0;j<K;++j){
            Out[i + j*K] = u[i] * v[j];
        }
    }
    return;
}

/* this works for matrices also if you just set the length "L" to be the right value,
 * often K*K, instead of just K in the case of vectors
 */
void componentVectorMultiplyInPlace(double * Out, double * u, double * v, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        Out[i] = u[i] * v[i];
    
    return;
}
