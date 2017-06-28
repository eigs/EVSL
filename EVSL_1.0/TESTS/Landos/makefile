.SUFFIXES : .c .o

include ../../makefile.in

FLAGS = -DUNIX -O3 -g -Wall

INCLUDES = -I../../INC  -I ../COMMON_GEN

EXTRA = ../COMMON_GEN/io.o ../COMMON_GEN/mmio.o ../COMMON_GEN/lapl.o

# Object files
OBJS_LD = LanDos.o exDOS.o  $(EXTRA)
OBJS_LDG = LanDosG.o exDOS.o  $(EXTRA)

LIB = -L../../ -levsl


LIB_EXT =  $(OMP) $(LIBLAPACK) $(LIB0)

# Rules
default: LanDos.ex LanDosG.ex
all: LanDos.ex LanDosG.ex

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
