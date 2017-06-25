.SUFFIXES : .c .o

include ../../makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../../INC -I../../EXTERNAL

ifeq ($(DIRECTSOL), SUITESPARSE)
  OBJ_DIRSOL = evsl_suitesparse.o
  include ../../EXTERNAL/makefile_suitesparse.in
  LIB_DIRSOL = $(LIB_UMF)
else
  OBJ_DIRSOL = evsl_cxsparse.o
  LIB_DIRSOL = -L../../EXTERNAL/CXSparse/Lib -lcxsparse
endif

# Object files
OBJS_RN = MMRLanN.o io.o mmio.o ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_RR = MMRLanR.o io.o mmio.o ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_PR = MMPLanR.o io.o mmio.o ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_PN = MMPLanN.o io.o mmio.o ../../EXTERNAL/$(OBJ_DIRSOL)

LIB = -L../../ -levsl
LIB_EXT = $(LIB_DIRSOL) $(OMP) $(LIBLAPACK) $(LIB0)

ALLEXE += MMRLanR.ex MMRLanN.ex MMPLanR.ex MMPLanN.ex

# Rules
default: MMRLanN.ex 

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe

MMRLanN.ex: $(OBJS_RN)
	$(LINK) -o MMRLanN.ex $(OBJS_RN) $(LIB) $(LIB_EXT)

MMRLanR.ex: $(OBJS_RR)
	$(LINK) -o MMRLanR.ex $(OBJS_RR) $(LIB) $(LIB_EXT)

MMPLanR.ex: $(OBJS_PR)
	$(LINK) -o MMPLanR.ex $(OBJS_PR) $(LIB) $(LIB_EXT)

MMPLanN.ex: $(OBJS_PN)
	$(LINK) -o MMPLanN.ex $(OBJS_PN) $(LIB) $(LIB_EXT)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
