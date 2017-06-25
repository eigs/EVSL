include ../makefile.in

INCLUDES = -I../INC -I.

ifeq ($(DIRECTSOL), SUITESPARSE)
  include makefile_suitesparse.in
  OBJ = evsl_suitesparse.o evsl_direct_f90.o
  INCLUDES += $(UMF_INC)
else
  OBJ = evsl_cxsparse.o evsl_direct_f90.o
  INCLUDES += -ICXSparse/Include
endif

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

all : $(OBJ)

clean:
	rm -f *.o
