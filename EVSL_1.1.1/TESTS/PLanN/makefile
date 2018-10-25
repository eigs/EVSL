include ../../makefile.in

INCLUDES = -I${EVSL_BUILD_DIR}/include -I../COMMON
EXTRA = ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o

# Object files
OBJS_L          = LapPLanN.o         $(EXTRA)
OBJS_MM         = MMPLanN.o          $(EXTRA)
OBJS_PN_MatFree = LapPLanN_MatFree.o $(EXTRA)

LIB = -L${EVSL_BUILD_DIR}/lib -levsl
LIB_EXT = $(LIBLAPACK) $(LIBBLAS) $(LIBS)

all: LapPLanN.ex LapPLanN_MatFree.ex MMPLanN.ex

# Exe
LapPLanN.ex: $(OBJS_L)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_L) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

LapPLanN_MatFree.ex: $(OBJS_PN_MatFree)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_PN_MatFree) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMPLanN.ex: $(OBJS_MM)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_MM) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

clean:
	rm -f *.o *.ex

distclean:
	rm -rf *.o *.ex OUT
