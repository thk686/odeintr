#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

R_CallMethodDef callMethods[]  = {{NULL, NULL, 0}};
R_CMethodDef cMethods[] = {{NULL, NULL, 0}};

void
R_init_myLib(DllInfo *info)
{
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}