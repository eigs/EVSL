.SUFFIXES : .c .o

include ../../makefile.in

FLAGS = -DUNIX -O3 -g -Wall

INCLUDES = -I../../INC -I../../EXTERNAL   -I ../COMMON_GEN

ifeq ($(DIRECTSOL), SUITESPARSE)
  OBJ_DIRSOL = evsl_suitesparse.o
  include ../../EXTERNAL/makefile_suitesparse.in
  LIB_DIRSOL = $(LIB_UMF)
else
  OBJ_DIRSOL = evsl_cxsparse.o
  LIB_DIRSOL = -L../../EXTERNAL/CXSparse/Lib -lcxsparse
  FLAGS += -DCXSPARSE=1
endif

EXTRA =  ../COMMON_GEN/io.o ../COMMON_GEN/mmio.o ../COMMON_GEN/lapl.o
# Object files
OBJS_KPM_RN = KPM_MMRLanN.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_KPM_RR = KPM_MMRLanR.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_KPM_PR = KPM_MMPLanR.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_KPM_PN = KPM_MMPLanN.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)

OBJS_Lap_RN = Lap_RLanN.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_Lap_RR = Lap_RLanR.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_Lap_PR = Lap_PLanR.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_Lap_PN = Lap_PLanN.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)

OBJS_Lan_RN = Lan_MMRLanN.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_Lan_RR = Lan_MMRLanR.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_Lan_PR = Lan_MMPLanR.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)
OBJS_Lan_PN = Lan_MMPLanN.o $(EXTRA) ../../EXTERNAL/$(OBJ_DIRSOL)


LIB = -L../../ -levsl
LIB_EXT = $(LIB_DIRSOL) $(OMP) $(LIBLAPACK) $(LIB0)

ALLEXE += KPM_MMRLanR.ex KPM_MMRLanN.ex KPM_MMPLanR.ex KPM_MMPLanN.ex
ALLEXE += Lan_MMRLanR.ex Lan_MMRLanN.ex Lan_MMPLanR.ex Lan_MMPLanN.ex
ALLEXE += Lap_RLanR.ex Lap_RLanN.ex Lap_PLanR.ex Lap_PLanN.ex

# Rules
default: KPM_MMRLanN.ex 

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe

KPM_MMRLanN.ex: $(OBJS_KPM_RN)
	$(LINK) -o KPM_MMRLanN.ex $(OBJS_KPM_RN) $(LIB) $(LIB_EXT)

KPM_MMRLanR.ex: $(OBJS_KPM_RR)
	$(LINK) -o KPM_MMRLanR.ex $(OBJS_KPM_RR) $(LIB) $(LIB_EXT)

KPM_MMPLanR.ex: $(OBJS_KPM_PR)
	$(LINK) -o KPM_MMPLanR.ex $(OBJS_KPM_PR) $(LIB) $(LIB_EXT)

KPM_MMPLanN.ex: $(OBJS_KPM_PN)
	$(LINK) -o KPM_MMPLanN.ex $(OBJS_KPM_PN) $(LIB) $(LIB_EXT)


Lap_RLanN.ex: $(OBJS_Lap_RN)
	$(LINK) -o Lap_RLanN.ex $(OBJS_Lap_RN) $(LIB) $(LIB_EXT)

Lap_RLanR.ex: $(OBJS_Lap_RR)
	$(LINK) -o Lap_RLanR.ex $(OBJS_Lap_RR) $(LIB) $(LIB_EXT)

Lap_PLanR.ex: $(OBJS_Lap_PR)
	$(LINK) -o Lap_PLanR.ex $(OBJS_Lap_PR) $(LIB) $(LIB_EXT)

Lap_PLanN.ex: $(OBJS_Lap_PN)
	$(LINK) -o Lap_PLanN.ex $(OBJS_Lap_PN) $(LIB) $(LIB_EXT)


Lan_MMRLanN.ex: $(OBJS_Lan_RN)
	$(LINK) -o Lan_MMRLanN.ex $(OBJS_Lan_RN) $(LIB) $(LIB_EXT)

Lan_MMRLanR.ex: $(OBJS_Lan_RR)
	$(LINK) -o Lan_MMRLanR.ex $(OBJS_Lan_RR) $(LIB) $(LIB_EXT)

Lan_MMPLanR.ex: $(OBJS_Lan_PR)
	$(LINK) -o Lan_MMPLanR.ex $(OBJS_Lan_PR) $(LIB) $(LIB_EXT)

Lan_MMPLanN.ex: $(OBJS_Lan_PN)
	$(LINK) -o Lan_MMPLanN.ex $(OBJS_Lan_PN) $(LIB) $(LIB_EXT)

all: $(ALLEXE)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
