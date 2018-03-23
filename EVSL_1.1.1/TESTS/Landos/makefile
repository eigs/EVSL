.SUFFIXES : .c .o

include ../../makefile.in

FLAGS ?= -DUNIX -O3 -g -Wall

GEN_INCLUDES = -I../../INC  -I ../COMMON_GEN
STD_INCLUDES = -I../../INC  -I ../COMMON

ifeq ($(INTEL), 1)
GEN_INCLUDES += -I${MKLROOT}/include
STD_INCLUDES += -I${MKLROOT}/include
endif

GEN_EXTRA = ../COMMON_GEN/io.o ../COMMON_GEN/mmio.o ../COMMON_GEN/lapl.o

STD_EXTRA = ../COMMON/io.o ../COMMON/mmio.o ../COMMON/lapl.o

# Object files
OBJS_LD = LanDos.o exDOS.o  $(STD_EXTRA)
OBJS_LDG = LanDosG.o exDOS.o  $(GEN_EXTRA)

LIB = -L../../ -levsl


LIB_EXT =  $(OMP) $(LIBLAPACK) $(LIB0)

# Rules
default: LanDos.ex LanDosG.ex
all: LanDos.ex LanDosG.ex

%G.o : %G.c
	$(CC) $(FLAGS) $(GEN_INCLUDES) -o $@ -x c -c $<

%.o : %.c
	$(CC) $(FLAGS) $(STD_INCLUDES) -o $@ -x c -c $<


# Exe
LanDos.ex: $(OBJS_LD)
	$(LINK) -o LanDos.ex $(OBJS_LD) $(LIB) $(LIB_EXT)

LanDosG.ex: $(OBJS_LDG)
	$(LINK) -o LanDosG.ex $(OBJS_LDG) $(LIB) $(LIB_EXT)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
