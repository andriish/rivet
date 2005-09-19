noinst_PROGRAMS = testrivet

myLDADD = $(top_builddir)/lib/libRivet.la 
myLDFLAGS = -export-dynamic

testrivet_SOURCES = testrivet.cc
testrivet_LDADD = $(myLDADD)
testrivet_LDFLAGS = $(myLDFLAGS)

check-local:
	./testrivet

