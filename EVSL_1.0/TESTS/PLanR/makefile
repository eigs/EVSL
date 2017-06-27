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
OBJS_L  =  LapPLanR.o $(EXTRA)
OBJS_MM  = MMPLanR.o $(EXTRA)
OBJS_MM_OMP = MMPLanR_omp.o $(EXTRA) 

LIB = -L../../ -levsl 
LIB_EXT = $(LIB_DIRSOL) $(OMP) $(LIBLAPACK) $(LIB0)

ALLEXE = LapPLanR.ex MMPLanR_omp.ex  MMPLanR.ex


# Rules
default: LapPLanR.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<
MMPLanR_omp.o : MMPLanR_omp.c
	$(CC) $(FLAGS) $(INCLUDES) $(OMP) -I ../../INC  -o $@ -c $<

# Exe
## polynomial
LapPLanR.ex: $(OBJS_L)
	$(LINK) -o LapPLanR.ex $(OBJS_L) $(LIB) $(LIBLAPACK) $(LIB0)


MMPLanR.ex: $(OBJS_MM)
	$(LINK) -o MMPLanR.ex $(OBJS_MM) $(LIB) $(LIB_EXT)

MMPLanR_omp.ex: $(OBJS_MM_OMP)
	$(LINK) -o MMPLanR_omp.ex $(OBJS_MM_OMP) $(LIB) $(LIB_EXT) $(OMP)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
