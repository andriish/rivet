// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/ExampleAnalysis.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {

  // Book histograms
  void ExampleAnalysis::init() {
    // Using histogram auto-booking is preferable if there are comparison datasets in HepData.
    // Since this is just a demo analysis, there is no associated paper!
    _histTot         = bookHistogram1D("TotalMult", "Total multiplicity", 100, -0.5, 99.5);
    _histChTot       = bookHistogram1D("TotalChMult", "Total charged multiplicity", 50, -1.0, 99.0);
    _histHadrTot     = bookHistogram1D("HadrTotalMult", "Total hadronic multiplicity", 100, -0.5, 99.5);
    _histHadrChTot   = bookHistogram1D("HadrTotalChMult", "Total hadronic charged multiplicity", 50, -1.0, 99.0);
    double edges[11] = { 0.5, 0.6, 0.7, 0.80, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0 };
    _histThrust      = bookHistogram1D("Thrust", "Thrust", vector<double>(edges, edges+11));
    _histMajor       = bookHistogram1D("Major", "Thrust major", 10, 0.0, 0.6);
    _histSphericity  = bookHistogram1D("Sphericity", "Sphericity", 10, 0.0, 0.8);
    _histAplanarity  = bookHistogram1D("Aplanarity", "Aplanarity", 10, 0.0, 0.3);
  }


  // Do the analysis
  void ExampleAnalysis::analyze(const Event& event) {
    Log log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    // Analyse and print some info
    const Multiplicity& cm = event.applyProjection(_cmultproj);
    const Multiplicity& cnm = event.applyProjection(_cnmultproj);
    log << Log::DEBUG << "Total multiplicity            = " << cnm.totalMultiplicity()  << endl;
    log << Log::DEBUG << "Total charged multiplicity    = " << cm.totalMultiplicity()   << endl;
    log << Log::DEBUG << "Hadron multiplicity           = " << cnm.hadronMultiplicity() << endl;
    log << Log::DEBUG << "Hadron charged multiplicity   = " << cm.hadronMultiplicity()  << endl;

    const Thrust& t = event.applyProjection(_thrustproj);
    log << Log::DEBUG << "Thrust = " << t.thrust() << endl;

    const Sphericity& s = event.applyProjection(_sphericityproj);
    log << Log::DEBUG << "Sphericity = " << s.sphericity() << endl;
    log << Log::DEBUG << "Aplanarity = " << s.aplanarity() << endl;

    // Fill histograms
    const double weight = event.weight();
    _histTot->fill(cnm.totalMultiplicity(), weight);
    _histChTot->fill(cm.totalMultiplicity(), weight);
    _histHadrTot->fill(cnm.hadronMultiplicity(), weight);
    _histHadrChTot->fill(cm.hadronMultiplicity(), weight);
    _histThrust->fill(t.thrust(), weight);
    _histMajor->fill(t.thrustMajor(), weight);
    _histSphericity->fill(s.sphericity(), weight);
    _histAplanarity->fill(s.aplanarity(), weight);

    // Finished
    log << Log::DEBUG << "Finished analyzing" << endl;
  }


  // Finalize
  void ExampleAnalysis::finalize() { 
    normalize(_histTot);
    normalize(_histChTot);
    normalize(_histHadrTot);
    normalize(_histHadrChTot);
    normalize(_histThrust);
    normalize(_histMajor);
    normalize(_histSphericity);
    normalize(_histAplanarity);
  }

}
