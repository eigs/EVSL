include ../makefile.in

DIRECT_DIRS = Interface
ifeq ($(DIRECTSOL), CXSPARSE)
   DIRECT_DIRS := CXSparse ${DIRECT_DIRS}
endif

all:
	@ \
	for i in ${DIRECT_DIRS}; \
	do \
	( cd $$i && $(MAKE) $@ ) || exit 1; \
	done

clean:
	@ \
	for i in ${DIRECT_DIRS}; \
	do \
	( cd $$i && $(MAKE) $@ ) || exit 1; \
	done

