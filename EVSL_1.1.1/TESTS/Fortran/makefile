include ../../makefile.in

# Object files
OBJ_DSOL = ${EVSL_SRC_DIR}/DIRECT/Interface/evsl_direct_f90.o ${DIRECTSOL_INTERFACE}
OBJS_PN = LapPLanN.o functns.o genmat.o
OBJS_PN_MatFree = LapPLanN_MatFree.o functns.o genmat.o
OBJS_PR = LapPLanR.o functns.o genmat.o
OBJS_RN = LapRLanN.o functns.o genmat.o $(OBJ_DSOL)
OBJS_RR = LapRLanR.o functns.o genmat.o $(OBJ_DSOL)
OBJ_CSR = csr_module.o

LIB = -L${EVSL_BUILD_DIR}/lib -levsl
LIB_EXT = $(LIBDIRECTSOL) $(LIBLAPACK_DIR) $(LIBLAPACK) $(LIBBLAS_DIR) $(LIBBLAS) $(LIBS_DIR) $(LIBS)
ifeq ($(FC),gfortran)
   LIB_EXT += -lgfortran
endif

all: LapPLanN.ex LapPLanR.ex LapPLanN_MatFree.ex LapRLanN.ex LapRLanR.ex

# a special rule
LapPLanN_MatFree.o: LapPLanN_MatFree.f90 $(OBJ_CSR)
	$(FC) $(FCFLAGS) $(INCLUDES) -o $@ -c LapPLanN_MatFree.f90

# Exe
LapPLanN.ex: $(OBJS_PN)
	$(LINK_FC) $(LDFLAGS) -o $@ $(OBJS_PN) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

LapPLanR.ex: $(OBJS_PR)
	$(LINK_FC) $(LDFLAGS) -o $@ $(OBJS_PR) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

LapPLanN_MatFree.ex: $(OBJS_PN_MatFree)
	$(LINK_FC) $(LDFLAGS) -o $@ $(OBJS_PN_MatFree) $(OBJ_CSR) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

LapRLanN.ex: $(OBJS_RN)
	$(LINK_FC) $(LDFLAGS) -o $@ $(OBJS_RN) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

LapRLanR.ex: $(OBJS_RR)
	$(LINK_FC) $(LDFLAGS) -o $@ $(OBJS_RR) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

clean:
	rm -f *.o *.mod *.ex

distclean:
	rm -rf *.o *.ex OUT
