include ../makefile.in
include makefile.in

OBJ = evsl_suitesparse.o evsl_suitesparse_f90.o

INCLUDES = -I../INC -I. $(UMF_INC)

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

all : $(OBJ)

clean:
	rm -f $(OBJ)
