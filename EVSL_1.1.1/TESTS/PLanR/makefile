.SUFFIXES : .c .o

include ../../makefile.in

FLAGS ?= -DUNIX -O3 -g -Wall 

INCLUDES = -I../../INC -I../../EXTERNAL   -I ../COMMON

EXTRA =  ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o
# Object files
OBJS_L  =  LapPLanR.o $(EXTRA)
OBJS_MM  = MMPLanR.o $(EXTRA)
OBJS_MM_OMP = MMPLanR_omp.o $(EXTRA) 

LIB = -L../../ -levsl 

ALLEXE = LapPLanR.ex MMPLanR_omp.ex  MMPLanR.ex


# Rules
default: LapPLanR.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c -x c $<
MMPLanR_omp.o : MMPLanR_omp.c
	$(CC) $(FLAGS) $(INCLUDES) $(OMP) -I ../../INC  -o $@ -c -x c $<

# Exe
## polynomial
LapPLanR.ex: $(OBJS_L)
	$(LINK) -o LapPLanR.ex $(OBJS_L) $(LIB) $(LIBLAPACK) $(LIB0)


MMPLanR.ex: $(OBJS_MM)
	$(LINK) -o MMPLanR.ex $(OBJS_MM) $(LIB) $(LIBLAPACK) $(LIB0)

MMPLanR_omp.ex: $(OBJS_MM_OMP)
	$(LINK) -o MMPLanR_omp.ex $(OBJS_MM_OMP) $(LIB) $(LIBLAPACK) $(LIB0) $(OMP)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
