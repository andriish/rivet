SUBDIRS = Examples

mySOURCES = AnalysisBase.cc RivetHandler.cc 

DOCFILES = AnalysisBase.h RivetHandler.h 

INCLUDEFILES = $(DOCFILES) \
               AnalysisBase.fh AnalysisBase.icc \
               RivetHandler.fh RivetHandler.icc

noinst_LTLIBRARIES = libRivetAnalysis.la
libRivetAnalysis_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
