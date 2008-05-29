check_PROGRAMS = testMatVec testBoost

testBoost_SOURCES = testBoost.cc
testBoost_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS)

testMatVec_SOURCES = testMatVec.cc
testMatVec_LDADD   = -lm -lgsl -lgslcblas
testMatVec_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS)

TESTS = testMatVec testBoost
