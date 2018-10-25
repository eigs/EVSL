include ../../makefile.in

INCLUDES = -I${EVSL_BUILD_DIR}/include  -I../COMMON
EXTRA = ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o

# Object files
OBJS_LD  = LanDos.o $(EXTRA)
OBJS_LDG = LanDosG.o $(EXTRA)

LIB = -L${EVSL_BUILD_DIR}/lib -levsl
LIB_EXT = $(LIBLAPACK_DIR) $(LIBLAPACK) $(LIBBLAS_DIR) $(LIBBLAS) $(LIBS)

all: LanDos.ex LanDosG.ex

# Exe
LanDos.ex: $(OBJS_LD)
	$(LINK_CC) $(LDFLAGS) -o LanDos.ex $(OBJS_LD) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

LanDosG.ex: $(OBJS_LDG)
	$(LINK_CC) $(LDFLAGS) -o LanDosG.ex $(OBJS_LDG) $(LIB) $(LIB_EXT)
	@mkdir -p OUT

clean:
	rm -f *.o *.ex

distclean:
	rm -rf *.o *.ex OUT
