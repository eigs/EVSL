.SUFFIXES : .c .o

include ../makefile.in

# Common flags
FLAGS = -DUNIX -O3 -g -Wall 

INCLUDES = -I../INC -ISRC 

# Object files
OBJS = 	evsl.o dumps.o cheblanNr.o cheblanTr.o chebpoly.o lanbounds.o misc_la.o \
	ratfilter.o ratlanNr.o ratlanTr.o spmat.o spslice.o vect.o chebsi.o     \

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
