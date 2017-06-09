.SUFFIXES : .c .o

include ../makefile.in

INCLUDES = -I../INC -ISRC 

# Object files
OBJS = 	evsl.o dumps.o cheblanNr.o cheblanTr.o chebpoly.o lanbounds.o \
        lanTrbounds.o misc_la.o ratfilter.o ratlanNr.o ratlanTr.o spmat.o \
        spslice.o vect.o chebsi.o landos.o spslice2.o simpson2.o lsPol1.o\
		dos_utils.o landosG2.o \
	../FORTRAN/evsl_f90.o

### OBJ_TIME = defined in makefile.in 
##mactime.o  
##timing.o

LIB = ../libevsl.a

# Rules
default: $(LIB)

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

$(LIB): $(OBJS) $(OBJ_TIME) 
	$(AR) $@ $(OBJS) $(OBJ_TIME)
	$(RANLIB) $@

clean:
	rm -f $(OBJS) $(OBJ_TIME) $(TIME_FILE) $(LIB)
