.SUFFIXES : .f .c .o

include ../makefile.in
include ../EXTERNAL/makefile.in

# Common flags
FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../INC -I../EXTERNAL $(UMF_INC)

# Object files for rational filter restart lanczos
OBJS_RR = MMRLanR.o io.o mmio.o ../EXTERNAL/evsl_suitesparse.o
# Object files for rational filter Non-restart lanczos
OBJS_RN = MMRLanN.o io.o mmio.o ../EXTERNAL/evsl_suitesparse.o
# Object files for polynomial filter restart lanczos
OBJS_PR = MMPLanR.o io.o mmio.o 
# Object files for polynomial filter restart lanczos (omp)
OBJS_PR_OMP = MMPLanR_omp.o io.o mmio.o 
# Object files for polynomial filter Non-restart lanczos
OBJS_PN = MMPLanN.o io.o mmio.o
# Object files for polynomial filter Non-restart lanczos
OBJS_PSI = MMPSI.o io.o mmio.o

LIB = -L.. -llancheb 

ALLEXE = MMPLanR.ex MMPLanN.ex MMPLanR_omp.ex \
	 MMRLanR.ex MMRLanN.ex MMPSI.ex

LIB_EXT = $(LIB_UMF) -fopenmp $(LIBLAPACK) $(LIB0)

# Rules
default: MMPLanN.ex 

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<
# Special rules for omp code
MMPLanR_omp.o : MMPLanR_omp.c
	$(CC) $(FLAGS) -I ../INC -fopenmp -o $@ -c $<

# Exe
MMPLanR.ex: $(OBJS_PR)
	$(LINK) -o MMPLanR.ex $(OBJS_PR) $(LIB) $(LIB_EXT)

MMPLanN.ex: $(OBJS_PN)
	$(LINK) -o MMPLanN.ex $(OBJS_PN) $(LIB) $(LIB_EXT)

MMPSI.ex: $(OBJS_PSI)
	$(LINK) -o MMPSI.ex $(OBJS_PSI) $(LIB) $(LIB_EXT)

MMPLanR_omp.ex: $(OBJS_PR_OMP)
	$(LINK) -o MMPLanR_omp.ex $(OBJS_PR_OMP) $(LIB) $(LIB_EXT) -fopenmp
	
MMRLanR.ex: $(OBJS_RR) 
	$(LINK) -o MMRLanR.ex $(OBJS_RR) $(LIB) $(LIB_EXT) 

MMRLanN.ex: $(OBJS_RN) 
	$(LINK) -o MMRLanN.ex $(OBJS_RN) $(LIB) $(LIB_EXT)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex 

cleanall:
	rm -f *.o *.ex OUT/*

