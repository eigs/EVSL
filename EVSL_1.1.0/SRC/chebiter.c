#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

/*!
 * @brief data needed for Chebyshev iterations
 *
 */
typedef struct _Chebiter_Data {
  /* eigenvalue bounds of the matrix */
  double lb, ub;
  /* polynomial degree */
  int deg;
  /* size */
  int n;
  /* matvec function and data */
  EVSLMatvec *mv;
  /* work space */
  double *w, *r, *p;
  /* residual norms */
  /*  @param: comp_res: 0: compute no residuals, 
   *                    1: compute the final residual
   *                    2: compute all residual history
   */
  int comp_res;
  double *res;
  /* stats */
  //size_t n_chebmv;
  //double t_chebmv;
} Chebiter_Data;

/** 
 * @brief Perform matrix-vector product y = A * x in Chebiter
 * 
 * */
static inline void EVSLChebMatvec(Chebiter_Data   *cheb_data, 
                                  double    *x, 
                                  double    *y) {

  //CHKERR(!cheb_data->mv);
     
  //double tms = pEVSL_Wtime();

  cheb_data->mv->func(x, y, cheb_data->mv->data);
  
  //double tme = pEVSL_Wtime();
  //cheb_data->t_chebmv += tme - tms;
  //cheb_data->n_chebmv ++;
}

/** @brief Setup Chebyshev iterations for a Parcsr matrix A for
 *         solving linear systems
 *  @param: tol, deg_in:
 *          when tol  > 0, degree will be found based on tol and 
 *                        deg_in serves as the max degree allowed;
 *          when tol <= 0, degree := degree_in
 * */
int EVSLChebIterSetup(double lmin, double lmax, csrMat *A,
                      double tol, int deg_in,
                      void **data) {
  Chebiter_Data *cheb;
  Malloc(cheb, 1, Chebiter_Data);
  int n = A->nrows, deg;
  double *w, *r, *p;
  Malloc(w, n, double);
  Malloc(r, n, double);
  Malloc(p, n, double);

  //cheb->n_chebmv = 0;
  //cheb->t_chebmv = 0.0;

  /* save the solver settings */
  cheb->lb  = lmin;
  cheb->ub  = lmax;
  cheb->n = n;
  Malloc(cheb->mv, 1, EVSLMatvec);
  cheb->mv->func = matvec_csr;
  cheb->mv->data = (void *) A;
  /* alloc work space */
  cheb->w = w;
  cheb->r = r;
  cheb->p = p;

  if (tol <= 0.0) {
    deg = deg_in;
  } else {
    EVSLChebIterFindDeg(tol, 10, cheb, &deg);
    deg = min(deg, deg_in);
  }
  deg = max(deg, 0);
  cheb->deg = deg;
  printf("ChebIter deg %d\n", deg);
  Calloc(cheb->res, deg+1, double);

  /* 0: do not compute any residual norm
   * 1: only save the final residual norm
   * 2: compute all residuals */
  cheb->comp_res = 0;

  *data = (void *) cheb;

  return 0;
}

/** @brief Solve function for Chebyshev iterations
 * Y. Saad, ``Iterative methods for sparse linear systems (2nd edition)'', 
 * Page 399
 */
void EVSLChebIterSolv(double *b, double *x, void *data) {
  int i;
  /* Cheb sol data */
  Chebiter_Data *Chebdata = (Chebiter_Data *) data;
  double theta, delta, alpha, beta, sigma, rho, rho1, tmp;
  double done=1.0, mdone = -1.0;
  int one = 1;
  double norm_r0, norm_r;
  double *res = Chebdata->res;

  int n = Chebdata->n;
  double *w = Chebdata->w;
  double *r = Chebdata->r;
  double *d = Chebdata->p;
  int deg = Chebdata->deg;
  int comp_res = Chebdata->comp_res;
  /* eig bounds */
  alpha = Chebdata->lb;
  beta  = Chebdata->ub;
  /* center and half width */
  theta = (beta + alpha) * 0.5;
  delta = (beta - alpha) * 0.5;
  sigma = theta / delta;
  rho   = 1.0 / sigma;
  /* use 0-initial guess, x_0 = 0, so r_0 = b */
  memset(x, 0, n*sizeof(double));
  memcpy(r, b, n*sizeof(double));
  /* d = 1/theta * r */
  memcpy(d, r, n*sizeof(double));
  tmp = 1.0 / theta;  DSCAL(&n, &tmp, d, &one);

  if (comp_res > 1) {
    norm_r0 = DNRM2(&n, r, &one);
    res[0] = norm_r0;
  }
  /* main iterations */
  for (i=0; i<deg; i++) {
    /* x = x + d */
    DAXPY(&n, &done, d, &one, x, &one);
    /* w = C * d */
    EVSLChebMatvec(Chebdata, d, w);
    /* r = r - w */
    DAXPY(&n, &mdone, w, &one, r, &one);
    /* rho1 = 1.0 / (2*sigma-rho) */
    rho1 = 1.0 / (2.0*sigma - rho);
    /* d = rho1*rho*d + 2*rho1/delta*r */
    tmp = rho1*rho;  DSCAL(&n, &tmp, d, &one);
    tmp = 2.0*rho1/delta;  DAXPY(&n, &tmp, r, &one, d, &one);
    /* update rho */
    rho = rho1;

    if (comp_res > 1) {
      norm_r = DNRM2(&n, r, &one);
      res[i+1] = norm_r;
    }
  }
  /* save the final residual norm */
  if (1 == comp_res) {
    res[0] = DNRM2(&n, r, &one);
  }
}

void EVSLChebIterFree(void *vdata) {
  Chebiter_Data *data = (Chebiter_Data *) vdata;
  free(data->w);
  free(data->r);
  free(data->p);
  free(data->mv);
  free(data->res);
  free(vdata);
}

/** @brief Find the Cheb poly degree for a given tolerance 
 */
void EVSLChebIterFindDeg(double tol, int nvec, void *data, int *degout) {
  Chebiter_Data *Chebdata = (Chebiter_Data *) data;
  double *b, *x;
  double theta, delta, alpha, beta, sigma, rho, rho1, tmp, tolr;
  double done=1.0, mdone = -1.0;
  int i, one=1, deg = 0;
  int n = Chebdata->n;
  double norm_r0, norm_r;
  double *w = Chebdata->w;
  double *r = Chebdata->r;
  double *d = Chebdata->p;
 
  Malloc(b, n, double);
  Malloc(x, n, double);

  /* average of nvec random vectors */
  randn_double(n, b);
  for (i=0; i<nvec-1; i++) {
    randn_double(n, x);
    DAXPY(&n, &done, x, &one, b, &one);
  }
  tmp = 1.0 / nvec;  DSCAL(&n, &tmp, b, &one);

  /* eig bounds */
  alpha = Chebdata->lb;
  beta  = Chebdata->ub;
  /* center and half width */
  theta = (beta + alpha) * 0.5;
  delta = (beta - alpha) * 0.5;
  sigma = theta / delta;
  rho   = 1.0 / sigma;
  /* use 0-initial guess, x_0 = 0, so r_0 = b */
  memset(x, 0, n*sizeof(double));
  memcpy(r, b, n*sizeof(double));
  /* d = 1/theta * r */
  memcpy(d, r, n*sizeof(double));
  tmp = 1.0 / theta;  DSCAL(&n, &tmp, d, &one);

  norm_r0 = DNRM2(&n, r, &one);
  tolr = tol * norm_r0;
  printf("ChebIter residual norm:\n");
  printf("  deg %3d: %e\n", deg, norm_r0);
  /* main iterations */
  while (1) {
    deg++;
    /* x = x + d */
    DAXPY(&n, &done, d, &one, x, &one);
    /* w = C * d */
    EVSLChebMatvec(Chebdata, d, w);
    /* r = r - w */
    DAXPY(&n, &mdone, w, &one, r, &one);
    /* rho1 = 1.0 / (2*sigma-rho) */
    rho1 = 1.0 / (2.0*sigma - rho);
    /* d = rho1*rho*d + 2*rho1/delta*r */
    tmp = rho1*rho;  DSCAL(&n, &tmp, d, &one);
    tmp = 2.0*rho1/delta;  DAXPY(&n, &tmp, r, &one, d, &one);
    /* update rho */
    rho = rho1;
    /* residual norm */
    norm_r = DNRM2(&n, r, &one);
    printf("  deg %3d: %e\n", deg, norm_r);
    if (norm_r < tolr) {
      break;
    }
  }

  *degout = deg;

  free(b);
  free(x);
}

