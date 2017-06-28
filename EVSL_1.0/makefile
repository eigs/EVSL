include makefile.in
all:
	( cd SRC                    ; make )
	( cd EXTERNAL               ; make )
ifeq ($(DIRECTSOL), CXSPARSE)
	( cd EXTERNAL/CXSparse      ; make )
endif
	( cd TESTS/Landos  ;  make all)
	( cd TESTS/Fortran ; make all)
	( cd TESTS/GEN ; make all)
	( cd TESTS/PLanR ; make all)
	( cd TESTS/RLanR ; make all)
	( cd TESTS/RLanN ; make all)
	( cd TESTS/PLanN ; make all)
	( cd TESTS/PSI ; make all)
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

docs:
	 ( doxygen Documentation/Doxyfile 2> Documentation/Doxygen-Errors.txt )
cleandocs:
	( cd Documentation; make  clean)

