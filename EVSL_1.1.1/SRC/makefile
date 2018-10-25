include ../makefile.in

SRC_DIRS=
ifeq ($(LIBBLAS), )
SRC_DIRS += blas
endif

ifeq ($(LIBLAPACK), )
SRC_DIRS += lapack
endif

INCLUDES = -I../INC
ifeq ($(INTEL_MKL), 1)
INCLUDES += -I${MKLROOT}/include
endif

# Object files
SRC = \
cheblanNr.c  chebsi.c     evsl.c         exDOS.c      landosG.c      ratfilter.c  simpson.c   spslice.c  vect.c \
cheblanTr.c  dos_utils.c  evsl_f90.c     lanbounds.c  lanTrbounds.c  ratlanNr.c   spmat.c     stats.c           \
chebpoly.c   dumps.c      evsl_memory.c  landos.c     misc_la.c      ratlanTr.c   spslice2.c  timing.c

OBJS = ${SRC:.c=.o}

all: $(OBJS)
	@ \
	for i in ${SRC_DIRS}; \
	do \
	   echo "Making $$i ..."; \
	   ( cd $$i && $(MAKE) $@) || exit 1; \
	   echo ""; \
	done


clean:
	@ \
	for i in ${SRC_DIRS}; \
	do \
	   echo "Cleaning $$i ..."; \
	   ( cd $$i && $(MAKE) $@ ) || exit 1; \
	   echo ""; \
	done
	rm -f *.o
