include ../makefile.in

INCLUDES = -I../INC -I.

ifeq ($(INTEL), 1)
  OBJ = evsl_pardiso.o 
  INCLUDES += -I${MKLROOT}/include
else
ifeq ($(DIRECTSOL), SUITESPARSE)
  include makefile_suitesparse.in
  OBJ = evsl_suitesparse.o
  INCLUDES += $(UMF_INC)
else
  OBJ = evsl_cxsparse.o
  INCLUDES += -ICXSparse/Include
endif
endif

OBJ += evsl_direct_f90.o

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c -x c $<

all : $(OBJ)

clean:
	rm -f *.o
