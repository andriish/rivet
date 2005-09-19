mySOURCES = Event.cc Projection.cc FinalStateProjection.cc

DOCFILES = Event.h Particle.h Projection.h FinalStateProjection.h

INCLUDEFILES = $(DOCFILES) Event.fh Event.icc Particle.fh \
               Projection.fh Projection.icc FinalStateProjection.icc

noinst_LTLIBRARIES = libRivetProjections.la
libRivetProjections_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
