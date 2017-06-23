include ../makefile.in
include makefile.in

OBJ += evsl_cxsparse.o
INCLUDES = -I../INC -I. -ICXSparse/Include

ifeq ($(USE_SUITESPARSE), 1)
OBJ += evsl_suitesparse.o evsl_suitesparse_f90.o
INCLUDES += $(UMF_INC)
endif

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

all : $(OBJ)

clean:
	rm -f $(OBJ)
