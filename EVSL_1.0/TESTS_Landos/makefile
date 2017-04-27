.SUFFIXES : .c .o

include ../makefile.in
include ../EXTERNAL/makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../INC -I../EXTERNAL $(UMF_INC)

# Object files
OBJS_LD = LanDos.o  ../EXTERNAL/evsl_suitesparse.o

LIB = -L../ -levsl

ALLEXE += LanDos.ex
LIB_EXT = $(LIB_UMF) -fopenmp $(LIBLAPACK) $(LIB0)

# Rules
default: LanDos.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe

LanDos.ex: $(OBJS_LD)
	$(LINK) -o LanDos.ex $(OBJS_LD) $(LIB) $(LIB_EXT)


all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
