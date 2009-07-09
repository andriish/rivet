// -*- C++ -*-
#include "Rivet/Analyses/STAR_2006_S6870392.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  STAR_2006_S6870392::STAR_2006_S6870392()
    : Analysis("STAR_2006_S6870392")
  {
    setBeams(PROTON, PROTON);

    /// @todo Use cross-section from generator
    //setNeedsCrossSection(true);

    //full final state
    FinalState fs(-2.0, 2.0);
    addProjection(fs, "FS");
    // R=0.4, pTmin=0, seed_threshold=0.5:
    addProjection(FastJets(fs, FastJets::CDFMIDPOINT, 0.4, 0.0, 0.5), "MidpointJets");
  } 


  // Book histograms
  void STAR_2006_S6870392::init() {
    _h_jet_pT_MB = 
      bookHistogram1D(1, 1, 1, "Inclusive jet cross-section, minbias trigger",
                      "jet $p_\\perp$/GeV","$1/(2\\pi) \\, d^2\\sigma/(d\\eta dp_\\perp)$ [pb/GeV]");
    _h_jet_pT_HT = 
      bookHistogram1D(2, 1, 1, "Inclusive jet cross-section, high tower trigger",
                      "jet $p_\\perp$/GeV","$1/(2\\pi) \\, d^2\\sigma/(d\\eta dp_\\perp)$ [pb/GeV]");
  }



  // Do the analysis 
  void STAR_2006_S6870392::analyze(const Event & event) {
    double weight = event.weight();

    // Skip if the event is empty
    const FinalState& fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because no final state found " << endl;
      vetoEvent;
    }

    // Find jets
    const FastJets& jetpro = applyProjection<FastJets>(event, "MidpointJets");
    const PseudoJets& jets = jetpro.pseudoJetsByPt();

    foreach (fastjet::PseudoJet jet, jets) {
      if (fabs(jets[0].eta()) < 0.8 && fabs(jets[0].eta()) > 0.2) {
        _h_jet_pT_MB->fill(jet.perp(), weight);
        _h_jet_pT_HT->fill(jet.perp(), weight);
      }
    }
  }



  // Finalize
  void STAR_2006_S6870392::finalize() {
    /// @todo Use the generator cross-section
    //_h_total_cross_section->fill(crossSection());
    normalize(_h_jet_pT_MB, 16603100);
    normalize(_h_jet_pT_HT, 1808234);
  }

}
