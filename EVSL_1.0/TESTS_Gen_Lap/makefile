.SUFFIXES : .c .o

include ../makefile.in
include ../EXTERNAL/makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../INC -I../EXTERNAL $(UMF_INC)

# Object files
OBJS_PR = LapPLanR.o io.o lapl.o ../EXTERNAL/evsl_suitesparse.o
OBJS_PN = LapPLanN.o io.o lapl.o ../EXTERNAL/evsl_suitesparse.o
OBJS_RN = LapRLanN.o io.o lapl.o ../EXTERNAL/evsl_suitesparse.o
OBJS_RR = LapRLanR.o io.o lapl.o ../EXTERNAL/evsl_suitesparse.o

LIB = -L../ -llancheb

ALLEXE += LapPLanR.ex LapPLanN.ex LapRLanN.ex LapRLanR.ex
LIB_EXT = $(LIB_UMF) -fopenmp $(LIBLAPACK) $(LIB0)

# Rules
default: LapPLanR.ex 

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe
LapPLanR.ex: $(OBJS_PR)
	$(LINK) -o LapPLanR.ex $(OBJS_PR) $(LIB) $(LIB_EXT)

LapPLanN.ex: $(OBJS_PN)
	$(LINK) -o LapPLanN.ex $(OBJS_PN) $(LIB) $(LIB_EXT)

LapRLanN.ex: $(OBJS_RN)
	$(LINK) -o LapRLanN.ex $(OBJS_RN) $(LIB) $(LIB_EXT)

LapRLanR.ex: $(OBJS_RR)
	$(LINK) -o LapRLanR.ex $(OBJS_RR) $(LIB) $(LIB_EXT)


all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
