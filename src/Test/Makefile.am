if WITH_TESTPROGS

noinst_PROGRAMS = testNewVector testBookFromAIDA

testBookFromAIDA_SOURCES  = testBookFromAIDA.cc
testBookFromAIDA_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)/include
testBookFromAIDA_LDADD    = -L$(top_builddir)/src -lRivet -L$(HEPMCLIBPATH) -l$(HEPMCLIBNAME)

testNewVector_SOURCES  = testNewVector.cc
testNewVector_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)/include/Rivet/Tools/D0RunIIcone

endif
