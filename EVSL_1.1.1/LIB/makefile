include ../makefile.in

all: libevsl${EVSL_LIB_SUFFIX}

BLAS_OBJS=../SRC/blas/*.o ../SRC/blas/f2c/*.o
LAPACK_OBJS=../SRC/lapack/*.o
OBJS=../SRC/*.o $(BLAS_OBJS) $(LAPACK_OBJS)

SONAME = libevsl-${EVSL_RELEASE_VERSION}${EVSL_LIB_SUFFIX}
SOLIBS = $(LIBLAPACKDIR) $(LIBLAPACK) $(LIBBLASDIR) $(LIBBLAS) ${LIBS}

libevsl.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

libevsl.so libevsl.dylib: ${OBJS}
	@echo  "Building $@ ... "
	${BUILD_CC_SHARED} -o ${SONAME} ${OBJS} ${SOLIBS} ${SHARED_SET_SONAME}${SONAME} ${SHARED_OPTIONS} ${LDFLAGS}
	ln -s ${SONAME} $@
	@ (cd .. && ln -fs LIB/${SONAME} $@) || exit 1

clean:
	rm -rf libevsl.*
