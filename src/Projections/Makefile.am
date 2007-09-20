LIBS = 
noinst_LTLIBRARIES = libRivetProjections.la
libRivetProjections_la_SOURCES  = \
  Beam.cc \
  ChargedFinalState.cc \
  ChargedLeptons.cc \
  CentralEtHCM.cc  \
  DISLepton.cc \
  DISKinematics.cc \
  D0ILConeJets.cc \
  FinalState.cc \
  FinalStateHCM.cc \
  Hemispheres.cc \
  KtJets.cc \
  Multiplicity.cc \
  ParisiTensor.cc \
  PVertex.cc \
  Sphericity.cc \
  Thrust.cc \
  TotalVisibleMomentum.cc \
  TrackJet.cc \
  VetoedFinalState.cc \
  WZandh.cc

libRivetProjections_la_CPPFLAGS = $(AM_CPPFLAGS)

if WITH_FASTJET
  libRivetProjections_la_SOURCES += FastJets.cc 
  libRivetProjections_la_CPPFLAGS += -DHAVE_FASTJET
endif
