.SUFFIXES : .c .o

include ../makefile.in

# Common flags
FLAGS = -DUNIX -O3 -g -Wall 

INC_UMF = -I$(SUITESPARSE_DIR)/UMFPACK/Include \
	  -I$(SUITESPARSE_DIR)/AMD/Include \
          -I$(SUITESPARSE_DIR)/SuiteSparse_config \
	  -I$(SUITESPARSE_DIR)/CHOLMOD/Include 
INCLUDES = -I../INC -ISRC 

# Object files
OBJS = 	vect.o cheblanTr.o cheblanNr.o ratlanTr.o ratlanNr.o ratfilter.o \
	misc_la.o lanbounds.o chebpoly.o spslice.o dumps.o \
	chebsi.o spmat.o evsl.o

ifneq ($(SUITESPARSE_DIR),)
  OBJS += suitesparse.o
  INCLUDES += $(INC_UMF)
  FLAGS += -DEVSL_WITH_SUITESPARSE
endif

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
