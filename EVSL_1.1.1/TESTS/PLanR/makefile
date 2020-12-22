include ../../makefile.in

INCLUDES = -I${EVSL_BUILD_DIR}/include -I ../COMMON
EXTRA =  ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o

ifeq ($(EVSL_WITH_CUDA), YES)
INCLUDES += -I${CUDA_HOME}/include
endif

# Object files
OBJS_L      = LapPLanR.o    $(EXTRA)
OBJS_MM     = MMPLanR.o     $(EXTRA)
OBJS_MM_OMP = MMPLanR_omp.o $(EXTRA)

LIB = -L${EVSL_BUILD_DIR}/lib -levsl
LIB_EXT = $(LIBLAPACK_DIR) $(LIBLAPACK) $(LIBBLAS_DIR) $(LIBBLAS) $(LIBS_DIR) $(LIBS)

all: LapPLanR.ex MMPLanR_omp.ex  MMPLanR.ex

# Exe
LapPLanR.ex: $(OBJS_L)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_L) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMPLanR.ex: $(OBJS_MM)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_MM) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMPLanR_omp.ex: $(OBJS_MM_OMP)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_MM_OMP) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

clean:
	rm -f *.o *.ex

distclean:
	rm -rf *.o *.ex OUT
