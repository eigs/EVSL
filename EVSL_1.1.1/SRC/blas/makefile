include ../../makefile.in

MUTE=@

# Object files
BLAS_SRC = \
daxpy.c   dcopy.c  dgemm.c  dger.c   dscal.c  dsymv.c  dsyr2k.c  dtrmv.c   lsame.c   zgemm.c  zscal.c  ztrsm.c \
dcabs1.c  ddot.c   dgemv.c  dnrm2.c  dswap.c  dsyr2.c  dtrmm.c   izamax.c  xerbla.c  zgeru.c  zswap.c

F2C_SRC = \
cabs.c  d_cnjg.c  d_imag.c  d_sign.c  i_nint.c  pow_di.c  s_cmp.c  s_copy.c  z_abs.c  z_div.c

BLAS_OBJS = ${BLAS_SRC:.c=.o}
F2C_OBJS = ${F2C_SRC:%.c=f2c/%.o}

all: $(BLAS_OBJS) ${F2C_OBJS}

INCLUDES = -I. -I../../INC

clean:
	rm -rf *.o f2c/*.o

distclean: clean
