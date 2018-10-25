include ../../makefile.in

INCLUDES = -I../../INC -I.. -I.

#ifeq ($(INTEL_MKL), 1)
#  OBJ = evsl_pardiso.o
#  INCLUDES += -I${MKLROOT}/include
#endif

ifeq ($(DIRECTSOL), SUITESPARSE)
INCLUDES += -I$(SUITESPARSE_DIR)/UMFPACK/Include -I$(SUITESPARSE_DIR)/AMD/Include -I$(SUITESPARSE_DIR)/SuiteSparse_config -I$(SUITESPARSE_DIR)/CHOLMOD/Include
endif

OBJ = ${DIRECTSOL_INTERFACE} evsl_direct_f90.o

all : $(OBJ)

clean:
	rm -f *.o
