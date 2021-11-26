#ifndef STRUCT_H
#define STRUCT_H

/**
 * @file struct.h
 * @brief structs used in evsl
 */

/**
 * @brief sparse matrix format: the coordinate (COO) format, 0-based
 *
 * ir, jc, vv : triples for all nonzeros (of size nnz)
 */
typedef struct _cooMat {
  int nrows,  /**< number of rows */
      ncols,  /**< number of columns */
      nnz,    /**< number of non-zeros */
      *ir,    /**< row indices of nonzero entries */
      *jc;    /**< column indices of a nonzero entries */
  double *vv; /**< values */
} cooMat;


/*!
 * @brief sparse matrix format: the compressed sparse row (CSR) format, 0-based
 *
 * 3-array variant: ia,ja,a, nnz == ia[nrows]
 */
typedef struct _csrMat {
  int owndata, /**< if owns (ia, ja, a) */
      /* on_dev, */  /**< NOT IN USE: RESERVED FOR LATER if (ia, ja, a) are device ptr */
      nrows,   /**< number of rows */
      ncols,   /**< number of columns */
      nnz,     /**< number of non-zeros */
      *ia,     /**< row pointers (of size nrows+1) */
      *ja;     /**< column indices (of size nnz) */
  double *a;   /**< numeric values (of size nnz) */
#ifdef EVSL_USING_CUDA_GPU
  cusparseMatDescr_t descr; /**< matrix descriptor for cusparse */
  char *dBuffer; /* GPU matvec may require buffer array */
#endif
} csrMat;

#ifdef EVSL_USING_CUDA_GPU

#if CUSPARSE_VERSION < 11000
#define EVSL_USING_CUSPARSE_HYB
#endif

#ifdef EVSL_USING_CUSPARSE_HYB
/*!
 * @brief cusparse HYB matrix format on GPU
 */
typedef struct _hybMat {
  int nrows,   /**< number of rows */
      ncols;   /**< number of columns */
  cusparseMatDescr_t descr;
  cusparseHybMat_t   hyb;
} hybMat;
#endif
#endif

/*!
 * @brief  parameters for polynomial filter
 *
 * default values are set by set_pol_def
 */
typedef struct _polparams {
  /** @name input to find_pol */
  /**@{*/
  int     max_deg;     /**< max allowed degree */
  int     min_deg ;    /**< min allowed degree */
  int     damping;     /**< 0 = no damping, 1 = Jackson, 2 = Lanczos */
  double  thresh_ext;  /**< threshold for accepting polynom. for end intervals */
  double  thresh_int;  /**< threshold for interior intervals */
  double  intvtol;     /**< cut-off point of middle interval */
  /**@}*/

  /** @name output from find_pol */
  /**@{*/
  int     type;       /**< type of the filter:
                           0: middle interval, 1: left interval, 2: right interval */
  double *mu;         /**< coefficients. allocation done by set_pol */
  double  cc;         /**< center of interval - used by chebAv */
  double  dd;         /**< half-width of interval - used by chebAv */
  double  gam;        /**< center of delta function used */
  double  bar;        /**< p(theta)>=bar indicates a wanted Ritz value */
  /**@}*/

  /** @name both input to and output from find_pol */
  /**@{*/
  int deg ;           /**< if deg == 0 before calling find_deg then
                       the polynomial degree is  computed
                       internally. Otherwise it is of degree deg.
                       [and  thresh_ext and thresh_int are not used]
                       default value=0, set by call to set_pol_def */
  /**@}*/
} polparams;

/**
 * @brief linear solver function prototype: [complex version]
 * which is used for solving system with A-SIGMA B
 * n  is the size  of the system,  br, bz are  the right-hand
 * side (real and  imaginary parts of complex vector),  xr, xz will
 * be the  solution (complex vector),  and "data" contains  all the
 * data  needed  by  the  solver.
 */
typedef void (*SolFuncC)(int n, double *br, double *bz, double *xr, double *xz, void *data);

/**
 * @brief function prototype for applying the solve B x = b
 */
typedef void (*SolFuncR)(double *b, double *x, void *data);

/**
 * @brief matvec function prototype
 */
typedef void (*MVFunc)(double *x, double *y, void *data);

/*!
 * @brief user-provided function and data for solving (A - SIGMA*B) x = b
 *
 */
typedef struct _EVSLASIGMABSol {
  SolFuncC func;       /**< function pointer */
  void *data;          /**< data */
} EVSLASIGMABSol;

/*!
 * @brief  parameters for rational filter
 *
 * default values are set by set_rat_def
 */
typedef struct _ratparams {
  /*  */
  int num;            /**< number of the poles */
  int pw;             /**< default multiplicity of each pole */
  int method;         /**< type of poles: 0: Cauchy pole, 1: Mid-point */
  double beta;        /**< LS approximation weight */
  double aa;          /**< left endpoint of the interval */
  double bb;          /**< right endpoint of the interval */
  double bar;         /**< rational function value at boundaries */
  /** internal data */
  int *mulp;          /**< multiplicity of each pole */
  int pow;            /**< total multiplicites of all poles */
  /** The following are output - i.e., set by find_ratf */
  EVSL_Complex *omega; /**< weights allocation done by find_ratf */
  EVSL_Complex *zk;    /**< locations of poles done by find_ratf */
  EVSLASIGMABSol *ASIGBsol; /**< function and data for A-&sigma B solve
                                 arrays of length ratparams.num */
} ratparams;


/*!
 * @brief user-provided Mat-Vec function and data for y = A * x or y = B * x
 *
 */
typedef struct _EVSLMatvec {
  MVFunc func;         /**< function pointer */
  void *data;          /**< data */
} EVSLMatvec;

/*!
 * @brief user-provided function and data for solving B x = b
 *
 */
typedef struct _EVSLBSol {
  SolFuncR func;       /**< function pointer */
  void *data;          /**< data */
} EVSLBSol;

/*!
 * @brief user-provided function for solving L^{T} x = b
 *
 */
typedef struct _EVSLLTSol {
  SolFuncR func;       /**< function pointer */
  void *data;          /**< data */
} EVSLLTSol;


/*!
 * @brief wrapper of all global variables in EVSL
 *
 */
typedef struct _evsldata {
  int n;                    /**< size of the problem [i.e., size(A), size(B)] */
  int ifGenEv;              /**< if it is a generalized eigenvalue problem */
  EVSLMatvec *Amv;          /**< matvec routine and the associated data for A */
  EVSLMatvec *Bmv;          /**< matvec routine and the associated data for B */
  EVSLBSol *Bsol;           /**< function and data for B solve */
  EVSLLTSol *LTsol;         /**< function and data for LT solve */
  double *ds;               /**< diagonal scaling matrix D,
                                 D^{-1}*A*D^{-1} = lambda * D^{-1}*B*D^{-1} */
#ifdef EVSL_USING_CUDA_GPU
  cublasHandle_t cublasH;
  cusparseHandle_t cusparseH;
  curandGenerator_t curandGen;
#endif
} evslData;

/*
 * Define a struct for L-S polynomial approximation to a matrix function
 */
typedef struct _BSolDataPol {
  int max_deg;    /**< max degree set by user */
  int deg;        /**< actual degree of the polynomial */
  double intv[2]; /**< spectrum bounds of the matrix */
  double cc, dd;  /**< center and half-width */
  double tol;     /**< approx. tolerance */
  double *mu;     /**< polyno. coeff */
  double *wk;     /**< working array for applying */
} BSolDataPol;


/* global variable: evslData */
extern evslData evsldata;

/*!
 * @brief timing and memory statistics of EVSL
 *
 */
typedef struct _evslstat {
  /* timing [level-0 funcs] */
  double t_setBsv;
  double t_setASigBsv;
  double t_iter;
  /* timing [level-1 funcs] */
  double t_polAv;
  double t_reorth;
  double t_siorth;
  double t_ritz;
  double t_eig;
  /* timing [level-2 funcs] */
  double t_mvA;
  double t_mvB;
  double t_svB;
  double t_svLT;
  double t_svASigB;
  double t_blas; //TODO
  double t_ratAv;
  size_t n_mvA;
  size_t n_mvB;
  size_t n_svB;
  size_t n_svLT;
  size_t n_svASigB;
  size_t n_polAv;
  size_t n_ratAv;
  /* memory TODO */
  size_t alloced;
  size_t alloced_total;
  size_t alloced_max;
} evslStat;

/* global variable: pevsl_stats */
extern evslStat evslstat;

#endif
