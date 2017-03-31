#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP odeintr_integrate_sys_const(SEXP derivsSEXP, SEXP recfunSEXP, SEXP initSEXP, SEXP durationSEXP, SEXP step_sizeSEXP, SEXP startSEXP, SEXP atolSEXP, SEXP rtolSEXP);
extern SEXP odeintr_integrate_sys_adapt(SEXP derivsSEXP, SEXP recfunSEXP, SEXP initSEXP, SEXP durationSEXP, SEXP step_sizeSEXP, SEXP startSEXP, SEXP atolSEXP, SEXP rtolSEXP);

R_CallMethodDef callMethods[]  = {
	{"odeintr_integrate_sys_const", (DL_FUNC) &odeintr_integrate_sys_const, 8},
	{"odeintr_integrate_sys_adapt", (DL_FUNC) &odeintr_integrate_sys_adapt, 8},
	{NULL, NULL, 0}};

void
R_init_myLib(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}