SUBDIRS = Examples

mySOURCES = AnalysisBase.cc RivetHandler.cc RivetInfo.cc

DOCFILES = AnalysisBase.h RivetHandler.h  RivetInfo.h

INCLUDEFILES = $(DOCFILES) AnalysisBase.fh AnalysisBase.icc \
               RivetHandler.fh RivetHandler.icc RivetInfo.icc

noinst_LTLIBRARIES = libRivetAnalysis.la
libRivetAnalysis_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
