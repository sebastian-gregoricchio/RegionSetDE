/* Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio. */

/* Registration of the C entry points, called from R through .Call() */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP C_readAlignmentRegions(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_cramTargets(SEXP);
extern SEXP C_cramHeadFlags(SEXP, SEXP);
extern SEXP C_cramIndexStats(SEXP);

static const R_CallMethodDef callMethods[] = {
    {"C_readAlignmentRegions", (DL_FUNC) &C_readAlignmentRegions, 9},
    {"C_cramTargets", (DL_FUNC) &C_cramTargets, 1},
    {"C_cramHeadFlags", (DL_FUNC) &C_cramHeadFlags, 2},
    {"C_cramIndexStats", (DL_FUNC) &C_cramIndexStats, 1},
    {NULL, NULL, 0}
};

void R_init_RegionSetDE(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
