include makefile.in
all:
	+ make -C SRC
	+ make -C EXTERNAL
ifeq ($(DIRECTSOL), CXSPARSE)
	+ make -C EXTERNAL/CXSparse
endif
	+ make all -C TESTS/Landos
	+ make all -C TESTS/GEN
	+ make all -C TESTS/PLanR
	+ make all -C TESTS/PLanN
	+ make all -C TESTS/RLanN
	+ make all -C TESTS/RLanR
	+ make all -C TESTS/PSI
	+ make all -C TESTS/Fortran

lib:
	+ make -C SRC
	+ make -C EXTERNAL
ifeq ($(DIRECTSOL), CXSPARSE)
	+ make -C EXTERNAL/CXSparse
endif

clean:
	( cd SRC                 ; make clean )
	( cd EXTERNAL            ; make clean )
ifeq ($(DIRECTSOL), CXSPARSE)
	( cd EXTERNAL/CXSparse   ; make clean )
endif
	( cd TESTS/Landos  ;  make cleanall)
	( cd TESTS/Fortran ; make clean)
	( cd TESTS/PLanR ; make cleanall)
	( cd TESTS/RLanR ; make cleanall)
	( cd TESTS/RLanN ; make cleanall)
	( cd TESTS/PLanN ; make cleanall)
	( cd TESTS/PSI ; make cleanall)
	( cd TESTS/GEN ; make cleanall)
	( cd TESTS/COMMON ; make clean)
	( cd TESTS/COMMON_GEN ; make clean)

docs:
	 ( doxygen Documentation/Doxyfile 2> Documentation/Doxygen-Errors.txt )
cleandocs:
	( cd Documentation; make  clean)

