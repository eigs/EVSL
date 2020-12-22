#include "cs.h"

void *_evsl_Malloc(size_t nbytes);
void *_evsl_Calloc(size_t count, size_t nbytes);
void *_evsl_Realloc(void *ptr, size_t nbytes);
void _evsl_Free(void *ptr);

/* wrapper for malloc */
void *cs_malloc (CS_INT n, size_t size)
{
    return (_evsl_Malloc(CS_MAX (n,1) * size)) ;
}

/* wrapper for calloc */
void *cs_calloc (CS_INT n, size_t size)
{
    return (_evsl_Calloc(CS_MAX (n,1), size)) ;
}

/* wrapper for free */
void *cs_free (void *p)
{
    if (p) _evsl_Free (p) ;  /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of cs_free */
}

/* wrapper for realloc */
void *cs_realloc (void *p, CS_INT n, size_t size, CS_INT *ok)
{
    void *pnew ;
    pnew = _evsl_Realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}

