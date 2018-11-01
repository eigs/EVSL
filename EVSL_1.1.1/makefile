MAKE_FILE_IN = makefile.in

ifneq ("$(wildcard $(MAKE_FILE_IN))","")
  include $(MAKE_FILE_IN)
else
  $(error $(MAKE_FILE_IN) does not exist. Run ./configure first)
endif

EVSL_DIRS = \
  SRC       \
  DIRECT    \
  LIB       \

EVSL_ALL_DIRS = ${EVSL_DIRS} TESTS
default: all

all:
	@ \
	for i in ${EVSL_DIRS}; \
	do \
	   echo "Making $$i ..."; \
	   ( cd $$i && $(MAKE) $@ ) || exit 1; \
	   echo ""; \
	done
	@ mkdir -p ${EVSL_BUILD_DIR}/include;
	@ mkdir -p ${EVSL_BUILD_DIR}/lib;
	cp -fR LIB/libevsl* ${EVSL_BUILD_DIR}/lib
	cp -fR ${EVSL_SRC_DIR}/INC/evsl.h ${EVSL_SRC_DIR}/INC/struct.h ${EVSL_SRC_DIR}/INC/EVSL_config.h ${EVSL_BUILD_DIR}/include

test: all
	@ \
	echo "Making test drivers ..."; \
	( cd TESTS; $(MAKE) clean; $(MAKE) )

install: all
	@ \
	echo "Installing EVSL in $(EVSL_INSTALL_DIR)"; \
	${EVSL_SRC_DIR}/CONFIG/mkinstalldirs ${EVSL_LIB_INSTALL} ${EVSL_INC_INSTALL}; \
	EVSL_PWD=`pwd`; \
	cd ${EVSL_BUILD_DIR}/lib; EVSL_FROMDIR=`pwd`; \
	cd $$EVSL_PWD; \
	cd ${EVSL_LIB_INSTALL}; EVSL_TODIR=`pwd`; \
	if [ "$$EVSL_FROMDIR" != "$$EVSL_TODIR" ]; \
	then \
	  cp -fR $$EVSL_FROMDIR/* $$EVSL_TODIR; \
	fi; \
	cd ${EVSL_BUILD_DIR}/include; EVSL_FROMDIR=`pwd`; \
	cd $$EVSL_PWD; \
	cd ${EVSL_INC_INSTALL}; EVSL_TODIR=`pwd`; \
	if [ "$$EVSL_FROMDIR" != "$$EVSL_TODIR" ]; \
	then \
	  cp -fR $$EVSL_FROMDIR/* $$EVSL_TODIR; \
	fi; \
	cd $$EVSL_PWD; \
	chmod -R a+rX,u+w,go-w ${EVSL_LIB_INSTALL}; \
	chmod -R a+rX,u+w,go-w ${EVSL_INC_INSTALL}; \
	echo

clean:
	@ \
	for i in ${EVSL_ALL_DIRS}; \
	   do \
	   echo "Cleaning $$i ..."; \
	   (cd $$i && $(MAKE) $@) || exit 1; \
	   echo ""; \
	done
	rm -rf ${EVSL_BUILD_DIR} libevsl.*

distclean: clean cleandocs
	rm -rf ./autom4te.cache ./config.log ./config.status makefile.in INC/EVSL_config.h

docs:
	( doxygen Documentation/Doxyfile 2> Documentation/Doxygen-Errors.txt )
cleandocs:
	( cd Documentation; make  clean)

