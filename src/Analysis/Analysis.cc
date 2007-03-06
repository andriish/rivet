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
#include "Rivet/Analysis/HepEx0112029.hh"
using namespace Rivet;

#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
using namespace AIDA;

#include <stdexcept>


//////////////////////////////////////////////////////////////


/// If you write a new analysis, add the ID for it here.
Analysis& Analysis::getAnalysis(const AnalysisName atype) {
  switch (atype) {
  case ANALYSIS_TEST:
    return *(new TestAnalysis());
  case ANALYSIS_HZ95108:
    return *(new HZ95108());
  case ANALYSIS_PL273B181:
    return *(new PL273B181());
  case ANALYSIS_HEPEX0112029:
    return *(new HepEx0112029());
  }
  throw runtime_error("Tried to get an analysis not known in the Rivet::AnalysisName enum.");
}


//////////////////////////////////////////////////////////////


Analysis::~Analysis() {}

void Analysis::init() {}


void Analysis::analyze(const Event &) {}


void Analysis::finalize() {}


RivetInfo Analysis::getInfo() const {
  return info;
}

//////////////////////////////////////////////////////////////


IAnalysisFactory& Analysis::analysisFactory() {
  return handler().analysisFactory();
}


ITree& Analysis::tree() {
  return handler().tree();
}


IHistogramFactory& Analysis::histogramFactory() {
  return handler().histogramFactory();
}


IHistogram1D* Analysis::bookHistogram1D(const std::string& name, const std::string& title, 
                                        const int nbins, const double lower, const double upper) {
  makeHistoDir();
  return histogramFactory().createHistogram1D(histoDir() + "/" + name, title, nbins, lower, upper);
}


void Analysis::makeHistoDir() {
  if (!_madeHistoDir) {
    if (!name().empty()) {
      tree().mkdir(histoDir());
    }
    _madeHistoDir = true;
  }
}
