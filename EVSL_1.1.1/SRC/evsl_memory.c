#include "internal_header.h"

void *_evsl_Malloc(size_t nbytes)
{
   void *ptr = malloc(nbytes);
   if (ptr == NULL)
   {
      fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", nbytes);
      CHKERR(1);
   }
   return ptr;
}

void *_evsl_Calloc(size_t count, size_t nbytes)
{
   void *ptr = calloc(count, nbytes);
   if (ptr == NULL)
   {
      fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", count * nbytes);
      CHKERR(1);
   }
   return ptr;
}

void *_evsl_Realloc(void *ptr, size_t nbytes)
{
   ptr = realloc(ptr, nbytes);
   if (ptr == NULL && nbytes > 0)
   {
      fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", nbytes);
      CHKERR(1);
   }
   return ptr;
}

void _evsl_Free(void *ptr)
{
   if (!ptr)
   {
      return;
   }

   free(ptr);
}

