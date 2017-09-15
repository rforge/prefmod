#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(calcs)(int *ncomp, int *notnaidx, int *iout);
extern void F77_NAME(calcs3)(int *ncomp, int *notnaidx, int *iout);

static const R_FortranMethodDef FortranEntries[] = {
  {"calcs",  (DL_FUNC) &F77_NAME(calcs),  3},
  {"calcs3", (DL_FUNC) &F77_NAME(calcs3), 3},
  {NULL, NULL, 0}
};

void R_init_prefmod(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
