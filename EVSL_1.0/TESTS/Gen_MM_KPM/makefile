.SUFFIXES : .c .o

include ../../makefile.in
include ../../EXTERNAL/makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../../INC -I../../EXTERNAL

# Object files
OBJS_RN = MMRLanN.o io.o mmio.o ../../EXTERNAL/evsl_cxsparse.o
OBJS_RR = MMRLanR.o io.o mmio.o ../../EXTERNAL/evsl_cxsparse.o
OBJS_PR = MMPLanR.o io.o mmio.o ../../EXTERNAL/evsl_cxsparse.o
OBJS_PN = MMPLanN.o io.o mmio.o ../../EXTERNAL/evsl_cxsparse.o

LIB = -L../../ -levsl

LIB_CXSPARSE=-L../../EXTERNAL/CXSparse/Lib -lcxsparse
ALLEXE += MMRLanR.ex MMRLanN.ex MMPLanR.ex MMPLanN.ex
LIB_EXT = $(LIB_CXSPARSE) $(OMP) $(LIBLAPACK) $(LIB0)

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
