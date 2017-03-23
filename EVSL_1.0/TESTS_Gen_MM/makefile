.SUFFIXES : .c .o

include ../makefile.in
include ../EXTERNAL/makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../INC -I../EXTERNAL $(UMF_INC)

# Object files
OBJS_RN = MMRLanN.o io.o mmio.o ../EXTERNAL/evsl_suitesparse.o
OBJS_RR = MMRLanR.o io.o mmio.o ../EXTERNAL/evsl_suitesparse.o

LIB = -L../ -llancheb

ALLEXE += MMRLanR.ex MMRLanN.ex 
LIB_EXT = $(LIB_UMF) -fopenmp $(LIBLAPACK) $(LIB0)

# Rules
default: MMRLanN.ex 

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe

MMRLanN.ex: $(OBJS_RN)
	$(LINK) -o MMRLanN.ex $(OBJS_RN) $(LIB) $(LIB_EXT)

MMRLanR.ex: $(OBJS_RR)
	$(LINK) -o MMRLanR.ex $(OBJS_RR) $(LIB) $(LIB_EXT)


all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
