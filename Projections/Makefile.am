mySOURCES = Event.cc Projection.cc BeamProjection.cc FinalStateProjection.cc \
            DISLepton.cc

DOCFILES = Event.h Particle.h Cmp.h Projection.h BeamProjection.h \
           FinalStateProjection.h DISLepton.h

INCLUDEFILES = $(DOCFILES) Event.fh Event.icc Particle.fh Cmp.fh Cmp.icc \
               Projection.fh Projection.icc BeamProjection.icc \
               FinalStateProjection.icc DISLepton.icc Particle.icc

noinst_LTLIBRARIES = libRivetProjections.la
libRivetProjections_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)
