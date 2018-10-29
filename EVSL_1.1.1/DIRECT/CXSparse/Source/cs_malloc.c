#include "cs.h"
#include "../../../INC/evsl.h"
#ifdef MATLAB_MEX_FILE
#define evsl_Malloc mxMalloc
#define evsl_Free mxFree
#define evsl_Realloc mxRealloc
#define evsl_Calloc mxCalloc
#endif

/* wrapper for malloc */
void *cs_malloc (CS_INT n, size_t size)
{
    return (evsl_Malloc (CS_MAX (n,1) * size, char)) ;
}

/* wrapper for calloc */
void *cs_calloc (CS_INT n, size_t size)
{
    return (evsl_Calloc (CS_MAX (n,1) * size, char)) ;
}

/* wrapper for free */
void *cs_free (void *p)
{
    if (p) evsl_Free (p) ;  /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of cs_free */
}

/* wrapper for realloc */
void *cs_realloc (void *p, CS_INT n, size_t size, CS_INT *ok)
{
    void *pnew ;
    pnew = evsl_Realloc (p, CS_MAX (n,1) * size, char) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}
