.SUFFIXES : .c .o

include ../../makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../../INC -I../../EXTERNAL   -I ../COMMON

ifeq ($(DIRECTSOL), SUITESPARSE)
  OBJ_DIRSOL = evsl_suitesparse.o
  include ../../EXTERNAL/makefile_suitesparse.in
  LIB_DIRSOL = $(LIB_UMF)
else
  OBJ_DIRSOL = evsl_cxsparse.o
  LIB_DIRSOL = -L../../EXTERNAL/CXSparse/Lib -lcxsparse
endif

EXTRA =  ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o
# Object files
OBJS_L  =  LapPLanN.o $(EXTRA)
OBJS_MM  = MMPLanN.o $(EXTRA)
OBJS_PN_MatFree = LapPLanN_MatFree.o $(EXTRA)

LIB = -L../../ -levsl 
LIB_EXT = $(LIB_DIRSOL) $(OMP) $(LIBLAPACK) $(LIB0)

ALLEXE = LapPLanN.ex LapPLanN_MatFree.ex  MMPLanN.ex


# Rules
default: LapPLanN.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe
## polynomial
LapPLanN.ex: $(OBJS_L)
	$(LINK) -o LapPLanN.ex $(OBJS_L) $(LIB) $(LIBLAPACK) $(LIB0)

LapPLanN_MatFree.ex: $(OBJS_PN_MatFree)
	$(LINK) -o LapPLanN_MatFree.ex $(OBJS_PN_MatFree) $(LIB) $(LIBLAPACK) $(LIB0)

MMPLanN.ex: $(OBJS_MM)
	$(LINK) -o MMPLanN.ex $(OBJS_MM) $(LIB) $(LIBLAPACK) $(LIB0)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
