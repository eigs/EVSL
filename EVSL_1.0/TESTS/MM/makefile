.SUFFIXES : .f .c .o

include ../../makefile.in

# Common flags
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

# Object files for rational filter restart lanczos
OBJS_RR = MMRLanR.o io.o mmio.o ../../EXTERNAL/$(OBJ_DIRSOL)
# Object files for rational filter Non-restart lanczos
OBJS_RN = MMRLanN.o io.o mmio.o ../../EXTERNAL/$(OBJ_DIRSOL)
# Object files for polynomial filter restart lanczos
OBJS_PR = MMPLanR.o io.o mmio.o 
# Object files for polynomial filter restart lanczos (omp)
OBJS_PR_OMP = MMPLanR_omp.o io.o mmio.o 
# Object files for polynomial filter Non-restart lanczos
OBJS_PN = MMPLanN.o io.o mmio.o
# Object files for polynomial filter Non-restart lanczos
OBJS_PSI = MMPSI.o io.o mmio.o

LIB = -L ../../ -levsl
LIB_EXT = $(LIB_DIRSOL) $(OMP) $(LIBLAPACK) $(LIB0)

ALLEXE = MMPLanR.ex MMPLanN.ex MMPLanR_omp.ex \
	 MMRLanR.ex MMRLanN.ex MMPSI.ex

# Rules
default: MMPLanN.ex 

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<
# Special rules for omp code
MMPLanR_omp.o : MMPLanR_omp.c
	$(CC) $(FLAGS) -I ../../INC $(OMP) -o $@ -c $<

# Exe
MMPLanR.ex: $(OBJS_PR)
	$(LINK) -o MMPLanR.ex $(OBJS_PR) $(LIB) $(LIB_EXT)

MMPLanN.ex: $(OBJS_PN)
	$(LINK) -o MMPLanN.ex $(OBJS_PN) $(LIB) $(LIB_EXT)

MMPSI.ex: $(OBJS_PSI)
	$(LINK) -o MMPSI.ex $(OBJS_PSI) $(LIB) $(LIB_EXT)

MMPLanR_omp.ex: $(OBJS_PR_OMP)
	$(LINK) -o MMPLanR_omp.ex $(OBJS_PR_OMP) $(LIB) $(LIB_EXT) $(OMP)

MMRLanR.ex: $(OBJS_RR) 
	$(LINK) -o MMRLanR.ex $(OBJS_RR) $(LIB) $(LIB_EXT) 

MMRLanN.ex: $(OBJS_RN) 
	$(LINK) -o MMRLanN.ex $(OBJS_RN) $(LIB) $(LIB_EXT)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex 

cleanall:
	rm -f *.o *.ex OUT/*

