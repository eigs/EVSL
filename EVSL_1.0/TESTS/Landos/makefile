.SUFFIXES : .c .o

include ../../makefile.in

FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../../INC 

# Object files
OBJS_LD = LanDos.o exDOS.o

LIB = -L../../ -levsl

## not needed for pol. filtering. 
LIB_UMF = 

LIB_EXT = $(LIB_UMF) -fopenmp $(LIBLAPACK) $(LIB0)

# Rules
default: LanDos.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe

LanDos.ex: $(OBJS_LD)
	$(LINK) -o LanDos.ex $(OBJS_LD) $(LIB) $(LIB_EXT)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
