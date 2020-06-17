#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "r2d2.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static const R_CallMethodDef callMethods[]  = {
  {"r2d2", (DL_FUNC) &r2d2, 5},
  {NULL, NULL, 0}
};

void R_init_C3PO(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

