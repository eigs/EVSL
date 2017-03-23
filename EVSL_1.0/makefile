ll:
	( cd SRC           ; make )
	( cd EXTERNAL      ; make )
	( cd TESTS_MM      ; mkdir -p OUT  ; make all)
	( cd TESTS_Lap     ; mkdir -p OUT  ; make all)
	( cd TESTS_Gen_Lap ; mkdir -p OUT  ; make all)
	( cd TESTS_Gen_MM  ; mkdir -p OUT  ; make all)
clean:
	( cd SRC        ; make clean)
	( cd EXTERNAL   ; make clean)
	( cd TESTS_Gen_Lap ; make cleanall)
	( cd TESTS_Gen_MM  ; make cleanall)	
	( cd TESTS_MM  ; make cleanall)
	( cd TESTS_Lap  ; make cleanall)

docs:
	 ( doxygen Documentation/Doxyfile 2> Documentation/Doxygen-Errors.txt )
cleandocs:
	( cd Documentation; make  clean)

