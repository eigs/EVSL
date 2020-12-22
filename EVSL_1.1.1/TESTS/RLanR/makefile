include ../../makefile.in

INCLUDES = -I${EVSL_BUILD_DIR}/include -I../../DIRECT/Interface -I ../COMMON
EXTRA =  ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o

ifeq ($(EVSL_WITH_CUDA), YES)
INCLUDES += -I${CUDA_HOME}/include
endif

# Object files
OBJS_L  = LapRLanR.o $(EXTRA) ${DIRECTSOL_INTERFACE}
OBJS_MM = MMRLanR.o  $(EXTRA) ${DIRECTSOL_INTERFACE}

LIB = -L${EVSL_BUILD_DIR}/lib -levsl
LIB_EXT = $(LIBDIRECTSOL) $(LIBLAPACK_DIR) $(LIBLAPACK) $(LIBBLAS_DIR) $(LIBBLAS) $(LIBS_DIR) $(LIBS)

all: LapRLanR.ex MMRLanR.ex

# Exe
LapRLanR.ex: $(OBJS_L)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_L) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMRLanR.ex: $(OBJS_MM)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_MM) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

clean:
	rm -f *.o *.ex

distclean:
	rm -rf *.o *.ex OUT
