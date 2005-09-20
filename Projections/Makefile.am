mySOURCES = Event.cc Projection.cc FinalStateProjection.cc DISLepton.cc

DOCFILES = Event.h Particle.h Projection.h FinalStateProjection.h DISLepton.h

INCLUDEFILES = $(DOCFILES) Event.fh Event.icc Particle.fh \
               Projection.fh Projection.icc FinalStateProjection.icc \
               DISLepton.icc Particle.icc

noinst_LTLIBRARIES = libRivetProjections.la
libRivetProjections_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
