#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_spacox_empirical_cgf_call(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"spacox_empirical_cgf", (DL_FUNC) &C_spacox_empirical_cgf_call, 4},
    {NULL, NULL, 0}
};

void attribute_visible R_init_SPACox(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
