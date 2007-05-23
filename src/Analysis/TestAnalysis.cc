// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/TestAnalysis.hh"
#include "Rivet/RivetAIDA.hh"
using namespace Rivet;

#include "AIDA/IHistogram1D.h"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;



/////////////////////////////////////////////////


// Book histograms
void TestAnalysis::init() {
  _histTot         = bookHistogram1D("TotalMult", "Total multiplicity", 100, -0.5, 999.5);
  _histChTot       = bookHistogram1D("TotalChMult", "Total charged multiplicity", 100, -0.5, 999.5);
  _histUnchTot     = bookHistogram1D("TotalUnchMult", "Total uncharged multiplicity", 100, -0.5, 999.5);
  _histHadrTot     = bookHistogram1D("HadrTotalMult", "Total hadronic multiplicity", 100, -0.5, 999.5);
  _histHadrChTot   = bookHistogram1D("HadrTotalChMult", "Total hadronic charged multiplicity", 100, -0.5, 999.5);
  _histHadrUnchTot = bookHistogram1D("HadrTotalUnchMult", "Total hadronic uncharged multiplicity", 100, -0.5, 999.5);
  _histThrust      = bookHistogram1D("Thrust", "Thrust", 100, 0.0, 1.0);
}


// Do the analysis
void TestAnalysis::analyze(const Event& event) {
  Log log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;

  // Analyse and print some info
  const Multiplicity& m = event.applyProjection(p_mult);
  log << Log::INFO << "Total multiplicity            = " << m.totalMultiplicity()           << endl;
  log << Log::INFO << "Total charged multiplicity    = " << m.totalChargedMultiplicity()    << endl;
  log << Log::INFO << "Total uncharged multiplicity  = " << m.totalUnchargedMultiplicity()  << endl;
  log << Log::DEBUG << "Hadron multiplicity           = " << m.hadronMultiplicity()          << endl;
  log << Log::DEBUG << "Hadron charged multiplicity   = " << m.hadronChargedMultiplicity()   << endl;
  log << Log::DEBUG << "Hadron uncharged multiplicity = " << m.hadronUnchargedMultiplicity() << endl;

  const Thrust& t = event.applyProjection(p_thrust);
  log << Log::INFO << "Thrust = " << t.thrust() << endl;

  // Fill histograms
  const double weight = event.weight();
  _histTot->fill(m.totalMultiplicity(), weight);
  _histChTot->fill(m.totalChargedMultiplicity(), weight);
  _histUnchTot->fill(m.totalUnchargedMultiplicity(), weight);
  _histHadrTot->fill(m.hadronMultiplicity(), weight);
  _histHadrChTot->fill(m.hadronChargedMultiplicity(), weight);
  _histHadrUnchTot->fill(m.hadronUnchargedMultiplicity(), weight);
  //
  _histThrust->fill(t.thrust(), weight);
  
  // Finished
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void TestAnalysis::finalize() { 
  normalize(_histTot);
  normalize(_histChTot);
  normalize(_histUnchTot);
  normalize(_histHadrTot);
  normalize(_histHadrChTot);
  normalize(_histHadrUnchTot);
  normalize(_histThrust);
}
