ll:
	( cd SRC    ; make )
	( cd TESTS_Gen  ; mkdir -p OUT  ; make all)
	( cd TESTS_Lap  ; mkdir -p OUT  ; make all)
clean:
	( cd SRC    ; make clean)
	( cd TESTS_Gen  ; make clean)
	( cd TESTS_Lap  ; make clean)

docs:
	 ( doxygen Documentation/Doxyfile 2> Documentation/Doxygen-Errors.txt )
cleandocs:
	( cd Documentation; make  clean)

