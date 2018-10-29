default: all

TEST_DIRS = \
 COMMON     \
 Fortran    \
 GEN        \
 Landos     \
 PLanN      \
 PLanR      \
 RLanN      \
 RLanR      \
 PSI

all:
	@ \
	for i in ${TEST_DIRS}; \
	do \
	   echo "Building $$i ..."; \
	   (cd $$i && $(MAKE) $@) || exit 1; \
	   echo ""; \
	done

clean:
	@ \
	for i in ${TEST_DIRS}; \
	do \
	   echo "Cleaning $$i ..."; \
	   (cd $$i && $(MAKE) $@) || exit 1; \
	   echo ""; \
	done

distclean:
	@ \
	for i in ${TEST_DIRS}; \
	do \
	   echo "Cleaning $$i ..."; \
	   (cd $$i && $(MAKE) $@) || exit 1; \
	   echo ""; \
	done

