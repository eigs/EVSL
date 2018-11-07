include ../../makefile.in

INCLUDES = -I${EVSL_BUILD_DIR}/include -I../../DIRECT/Interface -I../COMMON
EXTRA =  ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o

# Object files
OBJS_L  = LapRSI.o $(EXTRA) ${DIRECTSOL_INTERFACE}
OBJS_MM = MMRSI.o  $(EXTRA) ${DIRECTSOL_INTERFACE}

LIB = -L${EVSL_BUILD_DIR}/lib -levsl
LIB_EXT = $(LIBDIRECTSOL) $(LIBLAPACK_DIR) $(LIBLAPACK_DIR) $(LIBLAPACK) $(LIBBLAS_DIR) $(LIBBLAS) $(LIBS_DIR) $(LIBS)

all: LapRSI.ex MMRSI.ex

# Exe
LapRSI.ex: $(OBJS_L)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_L) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMRSI.ex: $(OBJS_MM)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_MM) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

clean:
	rm -f *.o *.ex

distclean:
	rm -rf *.o *.ex OUT
