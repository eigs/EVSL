.SUFFIXES : .c .o

include ../makefile.in

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

# Object files
OBJS_PN  = LapPLanN.o io.o
OBJS_PR  = LapPLanR.o io.o
OBJS_RN  = LapRLanN.o io.o
OBJS_RR  = LapRLanR.o io.o
OBJS_PSI = LapPSI.o io.o
OBJS_PN_MatFree = LapPLanN_MatFree.o io.o

LIB = -L../ -llancheb 

ALLEXE = LapPLanR.ex LapPLanN.ex LapPSI.ex LapPLanN_MatFree.ex
ifneq ($(SUITESPARSE_DIR),)
ALLEXE += LapRLanR.ex LapRLanN.ex
endif

# Rules
default: LapPLanN.ex

%.o : %.c
	$(CC) $(FLAGS)  -I ../INC -o $@ -c $<

# Exe
LapPLanR.ex: $(OBJS_PR) 
	$(LINK) -o LapPLanR.ex $(OBJS_PR) $(LIB) \
	$(LIBLAPACK) $(LIB0)

LapPLanN.ex: $(OBJS_PN) 
	$(LINK) -o LapPLanN.ex $(OBJS_PN) $(LIB) \
        $(LIBLAPACK) $(LIB0)

LapPSI.ex: $(OBJS_PSI) 
	$(LINK) -o LapPSI.ex $(OBJS_PSI) $(LIB) \
	$(LIBLAPACK) $(LIB0)

LapRLanR.ex: $(OBJS_RR) 
	$(LINK) -o LapRLanR.ex $(OBJS_RR) $(LIB) $(LIB_UMF) $(LIBLAPACK) \
	$(LIBBLAS) $(LIB0)

LapRLanN.ex: $(OBJS_RN) 
	$(LINK) -o LapRLanN.ex $(OBJS_RN) $(LIB) $(LIB_UMF) $(LIBLAPACK) \
	$(LIBBLAS) $(LIB0)

LapPLanN_MatFree.ex: $(OBJS_PN_MatFree)
	$(LINK) -o LapPLanN_MatFree.ex $(OBJS_PN_MatFree) $(LIB) \
        $(LIBLAPACK) $(LIB0)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
