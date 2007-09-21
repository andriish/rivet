// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/TestAnalysis.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {

  // Book histograms
  void TestAnalysis::init() {
    // Using histogram auto-booking is preferable if there are comparison datasets in HepData.
    // Since this is just a demo analysis, there is no associate paper!
    _histTot         = bookHistogram1D("TotalMult", "Total multiplicity", 100, -0.5, 99.5);
    _histChTot       = bookHistogram1D("TotalChMult", "Total charged multiplicity", 50, -1.0, 99.0);
    _histUnchTot     = bookHistogram1D("TotalUnchMult", "Total uncharged multiplicity", 100, -0.5, 99.5);
    _histHadrTot     = bookHistogram1D("HadrTotalMult", "Total hadronic multiplicity", 100, -0.5, 99.5);
    _histHadrChTot   = bookHistogram1D("HadrTotalChMult", "Total hadronic charged multiplicity", 50, -1.0, 99.0);
    _histHadrUnchTot = bookHistogram1D("HadrTotalUnchMult", "Total hadronic uncharged multiplicity", 100, -0.5, 99.5);
    //_histThrust      = bookHistogram1D("Thrust", "Thrust", 100, 0.5, 1.0);
    double edges[11] = { 0.0, 0.25, 0.5, 0.75, 0.80, 0.85, 0.9, 0.925, 0.95, 0.975, 1.0 };
    _histThrust     = bookHistogram1D("Thrust", "Thrust with different binning", vector<double>(edges, edges+11));
  }


  // Do the analysis
  void TestAnalysis::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    // Analyse and print some info
    const Multiplicity& m = event.applyProjection(_multproj);
    log << Log::DEBUG << "Total multiplicity            = " << m.totalMultiplicity()           << endl;
    log << Log::DEBUG << "Total charged multiplicity    = " << m.totalChargedMultiplicity()    << endl;
    log << Log::DEBUG << "Total uncharged multiplicity  = " << m.totalUnchargedMultiplicity()  << endl;
    log << Log::DEBUG << "Hadron multiplicity           = " << m.hadronMultiplicity()          << endl;
    log << Log::DEBUG << "Hadron charged multiplicity   = " << m.hadronChargedMultiplicity()   << endl;
    log << Log::DEBUG << "Hadron uncharged multiplicity = " << m.hadronUnchargedMultiplicity() << endl;

    const Thrust& t = event.applyProjection(_thrustproj);
    log << Log::DEBUG << "Thrust = " << t.thrust() << endl;

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

}
