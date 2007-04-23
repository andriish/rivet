if WITH_TESTPROGS

noinst_PROGRAMS = \
  testLogging testRivetHandler \
  testHepPDT testLWH testTinyXML \
  testD0RunIIcone testBookFromAIDA
MYLDADD  = -L$(top_builddir)/src -lRivet -L$(HEPMCLIBPATH) -l$(HEPMCLIBNAME)

testLogging_SOURCES  = testLogging.cc
testLogging_LDADD    = $(MYLDADD)

testRivetHandler_SOURCES  = testRivetHandler.cc
testRivetHandler_CPPFLAGS  = $(AM_CPPFLAGS) -I$(HEPMCINCPATH)
testRivetHandler_LDADD    = $(MYLDADD)

testHepPDT_SOURCES   = testHepPDT.cc
testHepPDT_CPPFLAGS  = $(AM_CPPFLAGS) -I$(HEPPIDINCPATH)
testHepPDT_LDADD     = $(MYLDADD)

testLWH_SOURCES  = testLWH.cc
testLWH_CPPFLAGS = $(AM_CPPFLAGS) -DLWH_USING_AIDA

testTinyXML_SOURCES  = testTinyXML.cc
testTinyXML_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)/include
testTinyXML_LDADD    = $(MYLDADD)

testBookFromAIDA_SOURCES  = testBookFromAIDA.cc
testBookFromAIDA_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)/include
testBookFromAIDA_LDADD    = $(MYLDADD)

testD0RunIIcone_SOURCES  = testD0RunIIcone.cc
testD0RunIIcone_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)/include/Rivet/Tools/D0RunIIcone


check-local:
	./testLogging
	./testRivetHandler
	./testHepPDT
	./testLWH && test -e test.aida.xml && rm -f test.aida.xml
#	./testTinyXML ../../data/4781.hepml.xml
	./testD0RunIIcone

endif