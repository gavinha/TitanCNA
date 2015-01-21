#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

/* Functions called by R code */

SEXP getPositionOverlapC(SEXP posn, SEXP start, SEXP stop);
SEXP fwd_backC_clonalCN(SEXP piGiZi, SEXP py, SEXP copyNumKey, SEXP zygosityKey, SEXP numClust, SEXP positions, SEXP zStrength, SEXP txnLen, SEXP useOutlier);
SEXP viterbiC_clonalCN(SEXP piGiZi, SEXP py, SEXP copyNumKey, SEXP zygosityKey, SEXP numClust, SEXP positions, SEXP zStrength, SEXP txnLen, SEXP useOutlier);
