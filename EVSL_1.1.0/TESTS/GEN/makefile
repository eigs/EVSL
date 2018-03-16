.SUFFIXES : .c .o

include ../../makefile.in

FLAGS ?= -DUNIX -O3 -g -Wall

INCLUDES = -I../../INC -I../../EXTERNAL   -I ../COMMON_GEN

ifeq ($(INTEL), 1)
  OBJ_DIRSOL = evsl_pardiso.o
else
ifeq ($(DIRECTSOL), SUITESPARSE)
  OBJ_DIRSOL = evsl_suitesparse.o
  include ../../EXTERNAL/makefile_suitesparse.in
  LIB_DIRSOL = $(LIB_UMF)
else
  OBJ_DIRSOL = evsl_cxsparse.o
  LIB_DIRSOL = -L../../EXTERNAL/CXSparse/Lib -lcxsparse
  FLAGS += -DCXSPARSE=1
endif
endif

EXTRA =  ../COMMON_GEN/io.o ../COMMON_GEN/mmio.o ../COMMON_GEN/lapl.o
# Object files
OBJS_S = SIMPLE.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_PN = MMPLanN.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_PR = MMPLanR.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_RN = MMRLanN.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_RR = MMRLanR.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)


LIB = -L../../ -levsl
LIB_EXT = $(LIB_DIRSOL) $(OMP) $(LIBLAPACK) $(LIB0)

ALLEXE = MMPLanR.ex MMPLanN.ex MMRLanR.ex MMRLanN.ex SIMPLE.ex

# Rules
default: SIMPLE.ex MMPLanN.ex MMRLanN.ex 

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -x c -c $<

##

MMRLanN.ex: $(OBJS_RN)
	$(LINK) -o MMRLanN.ex $(OBJS_RN) $(LIB) $(LIB_EXT)

MMRLanR.ex: $(OBJS_RR)
	$(LINK) -o MMRLanR.ex $(OBJS_RR) $(LIB) $(LIB_EXT)

MMPLanR.ex: $(OBJS_PR)
	$(LINK) -o MMPLanR.ex $(OBJS_PR) $(LIB) $(LIB_EXT)

MMPLanN.ex: $(OBJS_PN)
	$(LINK) -o MMPLanN.ex $(OBJS_PN) $(LIB) $(LIB_EXT)

SIMPLE.ex: $(OBJS_S)
	$(LINK) -o SIMPLE.ex $(OBJS_S) $(LIB) $(LIB_EXT)


all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
