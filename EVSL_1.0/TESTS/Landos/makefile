.SUFFIXES : .c .o

include ../../makefile.in

FLAGS = -DUNIX -O3 -g -Wall

INCLUDES = -I../../INC 

# Object files
OBJS_LD = LanDos.o exDOS.o 
OBJS_LDG = LanDosG.o exDOS.o io.o mmio.o 

LIB = -L../../ -levsl

## not needed for pol. filtering. 

LIB_EXT =  $(OMP) $(LIBLAPACK) $(LIB0)

# Rules
default: LanDos.ex LanDosG.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe
LanDos.ex: $(OBJS_LD)
	$(LINK) -o LanDos.ex $(OBJS_LD) $(LIB) $(LIB_EXT)

LanDosG.ex: $(OBJS_LDG)
	$(LINK) -o LanDosG.ex $(OBJS_LDG) $(LIB) $(LIB_EXT)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
