%module rivet

%{
  #define SWIG_FILE_WITH_INIT
  #include "Rivet/Analysis.hh"
  #include "Rivet/AnalysisHandler.hh"

  int run();
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_set.i"
%include "std_map.i"

namespace Rivet {
  class Projection;
}

%include "Rivet/ParticleName.hh"
%include "Rivet/HistoFormat.hh"
%include "Rivet/Cuts.fhh"
%include "Rivet/Cuts.hh"
%include "Rivet/Projection.fhh"
%include "Rivet/ProjectionApplier.hh"
%include "Rivet/Projection.hh"
%include "Rivet/Analysis.hh"
%include "Rivet/AnalysisHandler.hh"

int run();
