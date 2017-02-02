.SUFFIXES : .f .c .o

include ../makefile.in

# Common flags
FLAGS = -DUNIX -O3 -g -Wall 

LIB_METIS_5_1 = $(SUITESPARSE_DIR)/metis-5.1.0/build/Linux-x86_64/libmetis
LIB_UMF = -L$(SUITESPARSE_DIR)/UMFPACK/Lib -lumfpack \
          -L$(SUITESPARSE_DIR)/AMD/Lib/ -lamd \
          -L$(SUITESPARSE_DIR)/SuiteSparse_config -lsuitesparseconfig \
          -L$(SUITESPARSE_DIR)/CHOLMOD/Lib -lcholmod \
	  -L$(SUITESPARSE_DIR)/COLAMD/Lib -lcolamd \
          -L$(LIB_METIS_5_1) -lmetis \
          -L$(SUITESPARSE_DIR)/CCOLAMD/Lib -lccolamd \
	  -L$(SUITESPARSE_DIR)/CAMD/Lib -lcamd

INCLUDES = -I../INC  

UMF_INC = -I$(SUITESPARSE_DIR)/UMFPACK/Include \
	  -I$(SUITESPARSE_DIR)/AMD/Include \
          -I$(SUITESPARSE_DIR)/SuiteSparse_config \
          -I$(SUITESPARSE_DIR)/CHOLMOD/Include 

# Object files for rational filter restart lanczos
OBJS_RR = GenRLanR.o io.o mmio.o 
# Object files for rational filter Non-restart lanczos
OBJS_RN = GenRLanN.o io.o mmio.o 
# Object files for polynomial filter restart lanczos
OBJS_PR = GenPLanR.o io.o mmio.o 
# Object files for polynomial filter restart lanczos (omp)
OBJS_PR_OMP = GenPLanR_omp.o io.o mmio.o 
# Object files for polynomial filter Non-restart lanczos
OBJS_PN = GenPLanN.o io.o mmio.o
# Object files for subspace iteration
OBJS_PSI = GenPSI.o io.o mmio.o 

LIB = -L.. -llancheb 

ALLEXE = GenPLanR.ex GenPLanN.ex GenPSI.ex GenPLanR_omp.ex
ifneq ($(SUITESPARSE_DIR),)
ALLEXE += GenRLanR.ex GenRLanN.ex
endif

# Rules
default: GenPLanN.ex 

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<
# Special rules for omp code
GenPLanR_omp.o : GenPLanR_omp.c
	$(CC) $(FLAGS) -I ../INC -fopenmp -o $@ -c $<

# Exe
GenPLanR.ex: $(OBJS_PR)
	$(LINK) -o GenPLanR.ex $(OBJS_PR) $(LIB) $(LIBLAPACK) $(LIBBLAS) $(LIB0)

GenPLanN.ex: $(OBJS_PN)
	$(LINK) -o GenPLanN.ex $(OBJS_PN) $(LIB) $(LIBLAPACK) $(LIBBLAS) $(LIB0)

GenPSI.ex: $(OBJS_PSI)
	$(LINK) -o GenPSI.ex $(OBJS_PSI) $(LIB) $(LIBLAPACK) $(LIBBLAS) $(LIB0)

GenPLanR_omp.ex: $(OBJS_PR_OMP)
	$(LINK) -o GenPLanR_omp.ex $(OBJS_PR_OMP) $(LIB) $(LIBLAPACK) $(LIBBLAS) $(LIB0) -fopenmp
	
GenRLanR.ex: $(OBJS_RR) 
	$(LINK) -o GenRLanR.ex $(OBJS_RR) $(LIB) $(LIB_UMF) $(LIBLAPACK) $(LIBBLAS) $(LIB0)

GenRLanN.ex: $(OBJS_RN) 
	$(LINK) -o GenRLanN.ex $(OBJS_RN) $(LIB) $(LIB_UMF) $(LIBLAPACK) $(LIBBLAS) $(LIB0)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex 

cleanall:
	rm -f *.o *.ex OUT/*

