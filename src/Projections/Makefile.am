noinst_LTLIBRARIES = libRivetProjections.la
libRivetProjections_la_SOURCES = \
  Beam.cc \
  ChargedFinalState.cc \
  ChargedLeptons.cc \
  CentralEtHCM.cc \
  DISKinematics.cc \
  DISLepton.cc \
  D0ILConeJets.cc \
  FastJets.cc \
  FinalState.cc \
  FinalStateHCM.cc \
  HadronicFinalState.cc \
  Hemispheres.cc \
  JetShape.cc \
  LossyFinalState.cc \
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

if ENABLE_JADE
  libRivetProjections_la_CXXFLAGS = $(AM_CXXFLAGS) -D__HAVE_JADE
endif
