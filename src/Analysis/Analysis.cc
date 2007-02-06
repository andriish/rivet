// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Analysis abstract base class.
//

#include "Rivet/RivetHandler.hh"
#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Analysis/TestAnalysis.hh"
#include "Rivet/Analysis/HZ95108.hh"
#include "Rivet/Analysis/PL273B181.hh"
#include <stdexcept>

using namespace Rivet;


Analysis::~Analysis() {}


//////////////////////////////////////////////////////////////


Analysis* Analysis::getAnalysis(const AnalysisName atype) {
  switch (atype) {
  case ANALYSIS_TEST:
    return new TestAnalysis();
  case ANALYSIS_HZ95108:
    return new HZ95108();
  case ANALYSIS_PL273B181:
    return new PL273B181();
  }
  throw runtime_error("Tried to get an analysis not known in the Rivet::AnalysisName enum.");
}


//////////////////////////////////////////////////////////////


void Analysis::init() {}


void Analysis::analyze(const Event &) {}


void Analysis::finalize() {}


RivetInfo Analysis::getInfo() const {
  return info;
}

//////////////////////////////////////////////////////////////


AIDA::IAnalysisFactory* Analysis::analysisFactory() {
  return handler().analysisFactory();
}


AIDA::ITree* Analysis::tree() {
  return handler().tree();
}


AIDA::IHistogramFactory* Analysis::histogramFactory() {
  return handler().histogramFactory();
}

