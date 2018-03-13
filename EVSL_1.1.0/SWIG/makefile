default: lib

lib: lib_python

lib_python:
	python setup.py build_ext --inplace

clean_python:
	rm -rf build __pycache__ _evsl.* evsl.py evsl_wrap.c

clean: clean_python
