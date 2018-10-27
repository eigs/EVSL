/* INC/EVSL_config.h.  Generated from EVSL_config.h.in by configure.  */

/* Release name */
#define EVSL_RELEASE_NAME "EVSL"

/* Version number */
#define EVSL_RELEASE_VERSION "1.1.1"

/* Date of release */
#define EVSL_RELEASE_DATE "2018/03/23"

/* Time of release */
#define EVSL_RELEASE_TIME "00:00:00"

/* Bug reports */
#define EVSL_RELEASE_BUGS "saad@cs.umn.edu"

/* Define to 1 for Linux platforms */
#define EVSL_ARCH_LINUX 1

/* Define to 1 for Mac platforms */
/* #undef EVSL_ARCH_MAC */

/* Using internal BLAS routines in EVSL */
#define EVSL_USING_EVSL_BLAS 1

/* Using internal LAPACK routines in EVSL */
#define EVSL_USING_EVSL_LAPACK 1

/* Using Intel MKL for BLAS and LAPACK */
/* #undef EVSL_USING_INTEL_MKL */

/* Using long int version of Cholmod in SuiteSparse interface */
#define EVSL_CHOLMOD_USE_LONG 1

/* Enable OpenMP support */
/* #undef EVSL_USING_OPENMP */

/* EVSL basic data types */

/* EVSL data type for integer */
#define EVSL_Int int

/* EVSL data type for unsigned integer */
#define EVSL_Unsigned unsigned

/* EVSL data type for long integer */
#define EVSL_LongInt long int

/* EVSL data type for real */
#define EVSL_Real double

/* EVSL data type for complex */
#define EVSL_Complex double _Complex

/* Define to a macro mangling the given C identifier (in lower and upper
 * case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* EVSL Fortran 90 interface naming mangling */
#define EVSLFORT(name, NAME) FC_FUNC_(name ## _f90, NAME ## _f90)

