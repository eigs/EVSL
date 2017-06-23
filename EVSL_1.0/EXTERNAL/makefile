include ../makefile.in
include makefile.in

OBJ += evsl_cxsparse.o
INCLUDES = -I../INC -I. -ICXSparse/Include

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

all : $(OBJ)

clean:
	rm -f $(OBJ)
