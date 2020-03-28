#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _scBIRD_bird(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,SEXP);
extern SEXP _scBIRD_getLoci(SEXP);
extern SEXP _scBIRD_bird_select_loci(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_scBIRD_bird", (DL_FUNC) &_scBIRD_bird, 8},
    {"_scBIRD_getLoci", (DL_FUNC) &_scBIRD_getLoci, 1},
    {"_scBIRD_bird_select_loci", (DL_FUNC) &_scBIRD_bird_select_loci, 8},
    {NULL, NULL, 0}
};

void R_init_scBIRD(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
