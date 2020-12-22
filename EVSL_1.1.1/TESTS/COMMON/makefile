include ../../makefile.in

INCLUDES = -I../../INC  -I.

ifeq ($(EVSL_WITH_CUDA), YES)
INCLUDES += -I${CUDA_HOME}/include
endif

OBJ = io.o lapl.o mmio.o

all: ${OBJ}

clean:
	rm -f *.o

distclean: clean

