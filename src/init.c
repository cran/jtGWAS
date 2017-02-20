#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP jtGWAS_jtGWAS(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP jtGWAS_jtGWASmp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"jtGWAS_jtGWAS",   (DL_FUNC) &jtGWAS_jtGWAS,   5},
    {"jtGWAS_jtGWASmp", (DL_FUNC) &jtGWAS_jtGWASmp, 6},
    {NULL, NULL, 0}
};

void R_init_jtGWAS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

