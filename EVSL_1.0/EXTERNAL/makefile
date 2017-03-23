include ../makefile.in
include makefile.in

OBJ = evsl_suitesparse.o

INCLUDES = -I../INC -I. $(UMF_INC)

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

all : evsl_suitesparse.o

clean:
	rm -f $(OBJ)
