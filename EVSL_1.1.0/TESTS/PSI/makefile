.SUFFIXES : .c .o

include ../../makefile.in

FLAGS ?= -DUNIX -O3 -g -Wall 

INCLUDES = -I../../INC -I../../EXTERNAL   -I ../COMMON

EXTRA =  ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o
# Object files
OBJS_L  =  LapPSI.o $(EXTRA)
OBJS_MM  = MMPSI.o $(EXTRA)

LIB = -L../../ -levsl 
LIB_EXT = $(OMP) $(LIBLAPACK) $(LIB0)

ALLEXE = LapPSI.ex MMPSI.ex


# Rules
default: LapPSI.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -x c -c $<

# Exe
## polynomial
LapPSI.ex: $(OBJS_L)
	$(LINK) -o LapPSI.ex $(OBJS_L) $(LIB) $(LIB_EXT) 

MMPSI.ex: $(OBJS_MM)
	$(LINK) -o MMPSI.ex $(OBJS_MM) $(LIB) $(LIB_EXT) 

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
