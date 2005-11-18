noinst_PROGRAMS = testrivet testloop

myLDADD = $(top_builddir)/lib/libRivet.la -lCLHEP -lThePEG
myLDFLAGS = -export-dynamic \
            -L/usr/local/lib/ThePEG

testrivet_SOURCES = testrivet.cc
testrivet_LDADD = $(myLDADD)
testrivet_LDFLAGS = $(myLDFLAGS)

testloop_SOURCES = runEventLoop.cc
testloop_LDADD = $(myLDADD)
testloop_LDFLAGS = $(myLDFLAGS)

check-local:
	./testrivet
