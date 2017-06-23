.SUFFIXES : .c .o

include ../../makefile.in
include ../../EXTERNAL/makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../../INC -I../../EXTERNAL 

# Object files
OBJS_PN  = LapPLanN.o io.o lapl.o
OBJS_PR  = LapPLanR.o io.o lapl.o
OBJS_RN  = LapRLanN.o io.o lapl.o ../../EXTERNAL/evsl_cxsparse.o
OBJS_RR  = LapRLanR.o io.o lapl.o ../../EXTERNAL/evsl_cxsparse.o
OBJS_PSI = LapPSI.o io.o lapl.o
OBJS_PN_MatFree = LapPLanN_MatFree.o io.o lapl.o

LIB = -L../../ -levsl 

ALLEXE = LapPLanR.ex LapPLanN.ex LapPLanN_MatFree.ex LapRLanN.ex \
	 LapRLanR.ex LapPSI.ex

LIB_CXSPARSE=-L../../EXTERNAL/CXSparse/Lib -lcxsparse
LIB_EXT = $(LIB_CXSPARSE) $(OMP) $(LIBLAPACK) $(LIB0)

# Rules
default: LapPLanN.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe
## polynomial
LapPLanR.ex: $(OBJS_PR)
	$(LINK) -o LapPLanR.ex $(OBJS_PR) $(LIB) $(LIBLAPACK) $(LIB0)

LapPLanN.ex: $(OBJS_PN)
	$(LINK) -o LapPLanN.ex $(OBJS_PN) $(LIB) $(LIBLAPACK) $(LIB0)

LapPSI.ex: $(OBJS_PSI)
	$(LINK) -o LapPSI.ex $(OBJS_PSI) $(LIB) $(LIB_EXT)

LapPLanN_MatFree.ex: $(OBJS_PN_MatFree)
	$(LINK) -o LapPLanN_MatFree.ex $(OBJS_PN_MatFree) $(LIB) $(LIBLAPACK) $(LIB0)

## rational
LapRLanR.ex: $(OBJS_RR)
	$(LINK) -o LapRLanR.ex $(OBJS_RR) $(LIB) $(LIB_EXT)

LapRLanN.ex: $(OBJS_RN)
	$(LINK) -o LapRLanN.ex $(OBJS_RN) $(LIB) $(LIB_EXT)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
