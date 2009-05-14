// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/ExampleAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"

namespace Rivet {


  // Constructor
  ExampleAnalysis::ExampleAnalysis() {
    const FinalState cnfs;
    const ChargedFinalState cfs;
    addProjection(cnfs, "FS");
    addProjection(cfs, "CFS");
    addProjection(Multiplicity(cfs), "CMult");
    addProjection(Multiplicity(cnfs), "CNMult");
    addProjection(Thrust(cfs), "Thrust");
    addProjection(Sphericity(cfs), "Sphericity");
  }


  // Book histograms
  void ExampleAnalysis::init() {
    // Using histogram auto-booking is preferable if there are comparison datasets in HepData.
    // Since this is just a demo analysis, there is no associated paper!

    _histTot         = bookHistogram1D("TotalMult", "Total multiplicity", 
                                       "$N_\\text{tot}$", "$1/\\sigma \\, \\mathrm{d}{\\sigma}/\\mathrm{d}{N_\\text{tot}}$", 
                                       100, -0.5, 99.5);
    _histChTot       = bookHistogram1D("TotalChMult", "Total charged multiplicity", 
                                       "$N_\\text{ch}$", "$1/\\sigma \\, \\mathrm{d}{\\sigma}/\\mathrm{d}{N_\\text{ch}}$", 
                                       50, -1.0, 99.0);
    _histHadrTot     = bookHistogram1D("HadrTotalMult", "Total hadronic multiplicity", 
                                       "$N_\\text{H}$", "$1/\\sigma \\, \\mathrm{d}{\\sigma}/\\mathrm{d}{N_\\text{H}}$", 
                                       100, -0.5, 99.5);
    _histHadrChTot   = bookHistogram1D("HadrTotalChMult", "Total hadronic charged multiplicity", 
                                       "$N_\\text{Hch}$", "$1/\\sigma \\, \\mathrm{d}{\\sigma}/\\mathrm{d}{N_\\text{Hch}}$", 
                                       50, -1.0, 99.0);

    double edges[11] = { 0.5, 0.6, 0.7, 0.80, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0 };
    vector<double> vedges(edges, edges+11);
    _histThrust      = bookHistogram1D("Thrust", "Thrust", 
                                       "$T$", "$1/\\sigma \\, \\mathrm{d}{\\sigma}/\\mathrm{d}{T}$", vedges);
    _histMajor       = bookHistogram1D("Major", "Thrust major", 
                                       "$M$", "$1/\\sigma \\, \\mathrm{d}{\\sigma}/\\mathrm{d}{M}$", 10, 0.0, 0.6);
    _histSphericity  = bookHistogram1D("Sphericity", "Sphericity", 
                                       "$S$", "$1/\\sigma \\, \\mathrm{d}{\\sigma}/\\mathrm{d}{S}$", 10, 0.0, 0.8);
    _histAplanarity  = bookHistogram1D("Aplanarity", "Aplanarity",
                                       "$A$", "$1/\\sigma \\, \\mathrm{d}{\\sigma}/\\mathrm{d}{A}$", 10, 0.0, 0.3);
  }


  // Do the analysis
  void ExampleAnalysis::analyze(const Event& event) {
    // Analyse and print some info
    const Multiplicity& cm = applyProjection<Multiplicity>(event, "CMult");
    const Multiplicity& cnm = applyProjection<Multiplicity>(event, "CNMult");
    getLog() << Log::DEBUG << "Total multiplicity            = " << cnm.totalMultiplicity()  << endl;
    getLog() << Log::DEBUG << "Total charged multiplicity    = " << cm.totalMultiplicity()   << endl;
    getLog() << Log::DEBUG << "Hadron multiplicity           = " << cnm.hadronMultiplicity() << endl;
    getLog() << Log::DEBUG << "Hadron charged multiplicity   = " << cm.hadronMultiplicity()  << endl;

    const Thrust& t = applyProjection<Thrust>(event, "Thrust");
    getLog() << Log::DEBUG << "Thrust = " << t.thrust() << endl;

    const Sphericity& s = applyProjection<Sphericity>(event, "Sphericity");
    getLog() << Log::DEBUG << "Sphericity = " << s.sphericity() << endl;
    getLog() << Log::DEBUG << "Aplanarity = " << s.aplanarity() << endl;

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
