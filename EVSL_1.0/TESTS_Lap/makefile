.SUFFIXES : .c .o

include ../makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

LIB_METIS_5_1 = $(SUITESPARSE_DIR)/metis-5.1.0/build/Linux-x86_64/libmetis
LIB_UMF = -L$(SUITESPARSE_DIR)/UMFPACK/Lib -lumfpack \
          -L$(SUITESPARSE_DIR)/SuiteSparse_config -lsuitesparseconfig \
          -L$(SUITESPARSE_DIR)/CHOLMOD/Lib -lcholmod \
          -L$(SUITESPARSE_DIR)/AMD/Lib/ -lamd \
	  -L$(SUITESPARSE_DIR)/COLAMD/Lib -lcolamd \
          -L$(LIB_METIS_5_1) -lmetis \
          -L$(SUITESPARSE_DIR)/CCOLAMD/Lib -lccolamd \
	  -L$(SUITESPARSE_DIR)/CAMD/Lib -lcamd

INCLUDES = -I../INC 

UMF_INC = -I$(SUITESPARSE_DIR)/UMFPACK/Include \
	  -I$(SUITESPARSE_DIR)/AMD/Include \
          -I$(SUITESPARSE_DIR)/SuiteSparse_config \
          -I$(SUITESPARSE_DIR)/CHOLMOD/Include 

# Object files
OBJS_PN  = LapPLanN.o io.o lapl.o
OBJS_PR  = LapPLanR.o io.o lapl.o
OBJS_RN  = LapRLanN.o io.o lapl.o
OBJS_RR  = LapRLanR.o io.o lapl.o
OBJS_PSI = LapPSI.o io.o lapl.o
OBJS_PN_MatFree = LapPLanN_MatFree.o io.o lapl.o
OBJS_PR_GEN = LapPLanR_Gen.o io.o lapl.o
OBJS_RN_GEN = LapRLanN_Gen.o io.o lapl.o

LIB = -L../ -llancheb 

ALLEXE = LapPLanR.ex LapPLanN.ex LapPSI.ex LapPLanN_MatFree.ex
ifneq ($(SUITESPARSE_DIR),)
ALLEXE += LapRLanR.ex LapRLanN.ex LapPLanR_Gen.ex LapRLanN_Gen.ex
endif

ifneq ($(SUITESPARSE_DIR),)
LIB_EXT = $(LIB_UMF) -fopenmp
endif
LIB_EXT += $(LIBLAPACK) $(LIB0)

# Rules
default: LapPLanN.ex

%.o : %.c
	$(CC) $(FLAGS)  -I ../INC -o $@ -c $<

# Exe
## polynomial
LapPLanR.ex: $(OBJS_PR)
	$(LINK) -o LapPLanR.ex $(OBJS_PR) $(LIB) $(LIB_EXT) 	

LapPLanN.ex: $(OBJS_PN)
	$(LINK) -o LapPLanN.ex $(OBJS_PN) $(LIB) $(LIB_EXT)

LapPSI.ex: $(OBJS_PSI)
	$(LINK) -o LapPSI.ex $(OBJS_PSI) $(LIB) $(LIB_EXT)

LapPLanN_MatFree.ex: $(OBJS_PN_MatFree)
	$(LINK) -o LapPLanN_MatFree.ex $(OBJS_PN_MatFree) $(LIB) $(LIB_EXT)

## rational
LapRLanR.ex: $(OBJS_RR)
	$(LINK) -o LapRLanR.ex $(OBJS_RR) $(LIB) $(LIB_EXT)

LapRLanN.ex: $(OBJS_RN)
	$(LINK) -o LapRLanN.ex $(OBJS_RN) $(LIB) $(LIB_EXT)

## generalized e.v 
LapPLanR_Gen.ex: $(OBJS_PR_GEN)
	$(LINK) -o LapPLanR_Gen.ex $(OBJS_PR_GEN) $(LIB) $(LIB_EXT)

LapRLanN_Gen.ex: $(OBJS_RN_GEN)
	$(LINK) -o LapRLanN_Gen.ex $(OBJS_RN_GEN) $(LIB) $(LIB_EXT)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
