EXTRA_DIST = \
	__init__.py util.py \
    plotinfo.py spiresbib.py \
    core.pyx rivet.pxd \
    core.cpp


if WITH_CYTHON
core.cpp: core.pyx rivet.pxd
	cython core.pyx --cplus
else
core.cpp:
	@echo "Not (re)generating core.cpp since Cython is not installed"
endif


all-local: core.cpp
	@true

clean-local:
	rm -f core.cpp
