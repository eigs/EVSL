include ../../makefile.in

INCLUDES = -I${EVSL_BUILD_DIR}/include -I../../DIRECT/Interface -I../COMMON
EXTRA = ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o

ifeq ($(EVSL_WITH_CUDA), YES)
INCLUDES += -I${CUDA_HOME}/include
endif

# Object files
OBJS_S = MMsimple.o $(EXTRA) ${DIRECTSOL_INTERFACE}
OBJS_PN = MMPLanN.o $(EXTRA) ${DIRECTSOL_INTERFACE}
OBJS_PR = MMPLanR.o $(EXTRA) ${DIRECTSOL_INTERFACE}
OBJS_RN = MMRLanN.o $(EXTRA) ${DIRECTSOL_INTERFACE}
OBJS_RR = MMRLanR.o $(EXTRA) ${DIRECTSOL_INTERFACE}

LIB = -L${EVSL_BUILD_DIR}/lib -levsl
LIB_EXT = $(LIBDIRECTSOL) $(LIBLAPACK_DIR) $(LIBLAPACK) $(LIBBLAS_DIR) $(LIBBLAS) $(LIBS_DIR) $(LIBS)

all: MMPLanR.ex MMPLanN.ex MMRLanR.ex MMRLanN.ex MMsimple.ex

# Exe
MMRLanN.ex: $(OBJS_RN)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_RN) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMRLanR.ex: $(OBJS_RR)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_RR) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMPLanR.ex: $(OBJS_PR)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_PR) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMPLanN.ex: $(OBJS_PN)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_PN) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMsimple.ex: $(OBJS_S)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_S) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

clean:
	rm -f *.o *.ex

distclean:
	rm -rf *.o *.ex OUT
