include ../../makefile.in

INCLUDES = -I${EVSL_BUILD_DIR}/include -I../COMMON
EXTRA =  ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o

# Object files
OBJS_L  = LapPSI.o $(EXTRA)
OBJS_MM = MMPSI.o  $(EXTRA)

LIB = -L${EVSL_BUILD_DIR}/lib -levsl
LIB_EXT = $(LIBLAPACK) $(LIBBLAS) $(LIBS)

all: LapPSI.ex MMPSI.ex

# Exe
LapPSI.ex: $(OBJS_L)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_L) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

MMPSI.ex: $(OBJS_MM)
	$(LINK_CC) $(LDFLAGS) -o $@ $(OBJS_MM) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

clean:
	rm -f *.o *.ex

distclean:
	rm -rf *.o *.ex OUT
