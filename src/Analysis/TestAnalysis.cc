// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/TestAnalysis.hh"
using namespace Rivet;

#include "AIDA/IHistogram1D.h"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;



/////////////////////////////////////////////////


TestAnalysis::~TestAnalysis() {}


// Book histograms
void TestAnalysis::init() {
  histTot_         = bookHistogram1D("TotalMult", "Total multiplicity", 100, -0.5, 999.5);
  histChTot_       = bookHistogram1D("TotalChMult", "Total charged multiplicity", 100, -0.5, 999.5);
  histUnchTot_     = bookHistogram1D("TotalUnchMult", "Total uncharged multiplicity", 100, -0.5, 999.5);
  histHadrTot_     = bookHistogram1D("HadrTotalMult", "Total hadronic multiplicity", 100, -0.5, 999.5);
  histHadrChTot_   = bookHistogram1D("HadrTotalChMult", "Total hadronic charged multiplicity", 100, -0.5, 999.5);
  histHadrUnchTot_ = bookHistogram1D("HadrTotalUnchMult", "Total hadronic uncharged multiplicity", 100, -0.5, 999.5);
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

  // Fill histograms
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


// Provide info object
RivetInfo TestAnalysis::getInfo() const {
  return Analysis::getInfo() + fsproj.getInfo() + mult.getInfo();
}
