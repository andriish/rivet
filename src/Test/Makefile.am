if WITH_DEVTESTS

noinst_PROGRAMS = testMatVec testGSL

testMatVec_SOURCES = testMatVec.cc Vectors.hh Matrices.hh
testMatVec_LDADD   = -lm -lgsl -lgslcblas

testGSL_SOURCES = testGSL.cc
testGSL_LDADD   = -lgsl -lgslcblas -lm

endif
