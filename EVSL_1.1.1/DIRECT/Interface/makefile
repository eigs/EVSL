include ../../makefile.in

INCLUDES = -I../../INC -I.. -I.

ifeq ($(DIRECTSOL), MKL_PARDISO)
  INCLUDES += -I${MKLROOT}/include
endif

ifeq ($(DIRECTSOL), SUITESPARSE)
INCLUDES += -I$(SUITESPARSE_DIR)/UMFPACK/Include -I$(SUITESPARSE_DIR)/AMD/Include -I$(SUITESPARSE_DIR)/SuiteSparse_config -I$(SUITESPARSE_DIR)/CHOLMOD/Include
endif

ifeq ($(EVSL_WITH_CUDA), YES)
INCLUDES += -I${CUDA_HOME}/include
endif

all: ${DIRECTSOL_INTERFACE} evsl_direct_f90.o

clean:
	rm -f *.o
