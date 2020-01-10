#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(blackbox)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void F77_NAME(blackboxt)(int *, int *, int *, int *, double *, int *, int *, double *, double *, double *, double *, int *, int *, double *, int *);

extern void F77_NAME(mckalnew)(int *, int *, int *, int *, double *, int *, int *, double *, double *, double *, double *, double *, int *);  

static const R_FortranMethodDef FortranEntries[] = {
           {"blackbox", (DL_FUNC) &F77_NAME(blackbox), 15},
           {"blackboxt", (DL_FUNC) &F77_NAME(blackboxt), 15},
           {"mckalnew", (DL_FUNC) &F77_NAME(mckalnew), 13},
           {NULL, NULL, 0}
};

void R_init_basicspace(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
