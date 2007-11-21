noinst_LTLIBRARIES = libRivetProjections.la
libRivetProjections_la_SOURCES = \
  Beam.cc \
  ChargedFinalState.cc \
  ChargedLeptons.cc \
  CentralEtHCM.cc \
  DISKinematics.cc \
  DISLepton.cc \
  D0ILConeJets.cc \
  FinalState.cc \
  FinalStateHCM.cc \
  Hemispheres.cc \
  JetShape.cc \
  KtJets.cc \
  Multiplicity.cc \
  ParisiTensor.cc \
  PVertex.cc \
  Sphericity.cc \
  SVertex.cc \
  Thrust.cc \
  TotalVisibleMomentum.cc \
  TrackJet.cc \
  VetoedFinalState.cc \
  WZandh.cc 

libRivetProjections_la_CPPFLAGS = $(AM_CPPFLAGS)
## Make FastJet a definite requirement and ditch this conditional:
if WITH_FASTJET
  libRivetProjections_la_SOURCES += FastJets.cc 
  libRivetProjections_la_CPPFLAGS += -DHAVE_FASTJET
endif
