include ../../makefile.in
#include ../../EXTERNAL/makefile.in

FC = gfortran
FLAGS = -DUNIX -O3 -g

INCLUDES = -I../../INC -I../../EXTERNAL

ifeq ($(DIRECTSOL), SUITESPARSE)
  OBJ_DIRSOL = ../../EXTERNAL/evsl_suitesparse.o ../../EXTERNAL/evsl_direct_f90.o
  include ../../EXTERNAL/makefile_suitesparse.in
  LIB_DIRSOL = $(LIB_UMF)
else
  OBJ_DIRSOL = ../../EXTERNAL/evsl_cxsparse.o ../../EXTERNAL/evsl_direct_f90.o
  LIB_DIRSOL = -L../../EXTERNAL/CXSparse/Lib -lcxsparse
endif

# Object files
OBJS_PN = LapPLanN.o functns.o genmat.o
OBJS_PR = LapPLanR.o functns.o genmat.o
OBJS_PN_MatFree = LapPLanN_MatFree.o functns.o genmat.o csr_module.o
OBJS_RN = LapRLanN.o functns.o genmat.o $(OBJ_DIRSOL)
OBJS_RR = LapRLanR.o functns.o genmat.o $(OBJ_DIRSOL)

LIB = -L../../ -levsl $(LIBLAPACK)
LIB_EXT = $(LIB_DIRSOL) $(OMP) $(LIB0)

ALLEXE = LapPLanN.ex LapPLanR.ex LapPLanN_MatFree.ex LapRLanN.ex LapRLanR.ex


#Rules
default: LapPLanN.ex

#%.o : %.c
#	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

LapPLanN.o: LapPLanN.f90
	$(FC) $(FLAGS) -o $@ -c $<

LapPLanR.o: LapPLanR.f90
	$(FC) $(FLAGS) -o $@ -c $<

LapPLanN_MatFree.o: LapPLanN_MatFree.f90 csr_module.o csrmatrix.mod
	$(FC) $(FLAGS) -o $@ -c $<

LapRLanN.o: LapRLanN.f90
	$(FC) $(FLAGS) -o $@ -c $<

LapRLanR.o: LapRLanR.f90
	$(FC) $(FLAGS) -o $@ -c $<

functns.o: functns.f90
	$(FC) $(FLAGS) -o $@ -c $<

csr_module.o: csr_module.f90
	$(FC) $(FLAGS) -o $@ -c $<

genmat.o: genmat.f
	$(FC) $(FLAGS) -o $@ -c $<

# Exe
# Polynomial
LapPLanN.ex: $(OBJS_PN)
	$(LINK) -o LapPLanN.ex $(OBJS_PN) $(LIB) $(LIB0)

LapPLanR.ex: $(OBJS_PR)
	$(LINK) -o LapPLanR.ex $(OBJS_PN) $(LIB) $(LIB0)

LapPLanN_MatFree.ex: $(OBJS_PN_MatFree)
	$(LINK) -o LapPLanN_MatFree.ex $(OBJS_PN_MatFree) $(LIB) $(LIB0)

LapRLanN.ex: $(OBJS_RN)
	$(LINK) -o LapRLanN.ex $(OBJS_RN) $(LIB_EXT) $(LIB)

LapRLanR.ex: $(OBJS_RR)
	$(LINK) -o LapRLanR.ex $(OBJS_RR) $(LIB_EXT) $(LIB)

all: $(ALLEXE)

clean:
	rm -f *.o *.mod *.ex
