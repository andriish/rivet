// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/TestAnalysis.hh"
#include "HepPDT/ParticleID.hh"

using namespace Rivet;
using namespace HepMC;
using namespace AIDA;

TestAnalysis::~TestAnalysis() {}


// Book histograms
void TestAnalysis::init() {

  // Book histograms
  IHistogramFactory* hf = histogramFactory();
  tree()->mkdir("/Test");
  histTot_         = hf->createHistogram1D("/Test/TotalMult", 
                                           "Total multiplicity", 100, -0.5, 999.5);
  histChTot_       = hf->createHistogram1D("/TotalChMult", 
                                           "Total charged multiplicity", 51, -0.5, 50.5);
  histUnchTot_     = hf->createHistogram1D("/TotalUnchMult", 
                                           "Total uncharged multiplicity", 51, -0.5, 50.5);
  histHadrTot_     = hf->createHistogram1D("/HadrTotalMult", 
                                           "Total hadronic multiplicity", 51, -0.5, 50.5);
  histHadrChTot_   = hf->createHistogram1D("/HadrTotalChMult", 
                                           "Total hadronic charged multiplicity", 51, -0.5, 50.5);
  histHadrUnchTot_ = hf->createHistogram1D("/HadrTotalUnchMult", 
                                           "Total hadronic uncharged multiplicity", 51, -0.5, 50.5);
}


// Do the analysis
void TestAnalysis::analyze(const Event & event) {
  Logger& log = getLogger(ANALYSIS_TEST);
  log << LogPriority::DEBUG << "Starting analyzing" << endlog;

  // Analyse and print some info
  const Multiplicity& m = event.applyProjection(mult);
  log << LogPriority::INFO << "Total multiplicity            = " << m.totalMultiplicity()           << endlog;
  log << LogPriority::INFO << "Total charged multiplicity    = " << m.totalChargedMultiplicity()    << endlog;
  log << LogPriority::INFO << "Total uncharged multiplicity  = " << m.totalUnchargedMultiplicity()  << endlog;
  log << LogPriority::INFO << "Hadron multiplicity           = " << m.hadronMultiplicity()          << endlog;
  log << LogPriority::INFO << "Hadron charged multiplicity   = " << m.hadronChargedMultiplicity()   << endlog;
  log << LogPriority::INFO << "Hadron uncharged multiplicity = " << m.hadronUnchargedMultiplicity() << endlog;

  // Fill histograms here
  histTot_->fill(m.totalMultiplicity(), 1.0);
  histChTot_->fill(m.totalChargedMultiplicity(), 1.0);
  histUnchTot_->fill(m.totalUnchargedMultiplicity(), 1.0);
  histHadrTot_->fill(m.hadronMultiplicity(), 1.0);
  histHadrChTot_->fill(m.hadronChargedMultiplicity(), 1.0);
  histHadrUnchTot_->fill(m.hadronUnchargedMultiplicity(), 1.0);
  
  // Finished...
  log << LogPriority::DEBUG << "Finished analyzing" << endlog;
}


// Finalize
void TestAnalysis::finalize() { }


RivetInfo TestAnalysis::getInfo() const {
  return Analysis::getInfo() + mult.getInfo();
}
