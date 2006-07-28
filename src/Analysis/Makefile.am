noinst_LTLIBRARIES = libRivetAnalysis.la libRivetAnalysisExamples.la

mySOURCES1    = AnalysisBase.cc RivetHandler.cc RivetInfo.cc
DOCFILES1     = AnalysisBase.h RivetHandler.h  RivetInfo.h
INCLUDEFILES1 = $(DOCFILES1) AnalysisBase.fh AnalysisBase.icc \
               RivetHandler.fh RivetHandler.icc RivetInfo.icc
libRivetAnalysis_la_SOURCES = $(mySOURCES1) $(INCLUDEFILES1)

mySOURCES2    = Examples/TestMultiplicity.cc Examples/HZ95108.cc
DOCFILES2     = Examples/TestMultiplicity.h Examples/HZ95108.h
INCLUDEFILES2 = $(DOCFILES2) Examples/TestMultiplicity.icc Examples/HZ95108.icc
libRivetAnalysisExamples_la_SOURCES = $(mySOURCES2) $(INCLUDEFILES2)
