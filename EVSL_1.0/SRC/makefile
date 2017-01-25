.SUFFIXES : .c .o

include ../makefile.in

# Common flags
FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../INC -ISRC

# Object files
OBJS = 	vect.o cheblanTr.o cheblanNr.o ratlanTr.o ratlanNr.o ratfilter.o misc_la.o lanbounds.o \
	chebpoly.o spslice.o exeiglap3.o dumps.o chebsi.o spmat.o

### OBJ_TIME = defined in makefile.in 
##mactime.o  
##timing.o

LIB = ../liblancheb.a

# Rules
default: $(LIB)

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<
$(LIB): $(OBJS) $(OBJ_TIME) 
	$(AR) $@ $(OBJS) $(OBJ_TIME)
	$(RANLIB) $@

clean:
	rm -f $(OBJS) $(OBJ_TIME) $(TIME_FILE) $(LIB)
