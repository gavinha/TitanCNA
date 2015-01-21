#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "titan.h"
         
//Register "forward_backward" and "viterbi"
static const R_CallMethodDef callMethods[]  = {
  {"getPositionOverlapC", (DL_FUNC) &getPositionOverlapC, 3},
  {"fwd_backC_clonalCN", (DL_FUNC) &fwd_backC_clonalCN, 9},
  {"viterbiC_clonalCN", (DL_FUNC) &viterbiC_clonalCN, 9},
  {NULL, NULL, 0}
};

void R_init_TitanCNA(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
