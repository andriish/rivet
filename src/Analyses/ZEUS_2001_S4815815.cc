// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/ZEUS_2001_S4815815.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  // Constructor
  ZEUS_2001_S4815815::ZEUS_2001_S4815815() { 
      setBeams(POSITRON, PROTON);
      FinalState fs;
      addProjection(fs, "FS");
      /// @todo This is the *wrong* jet def: correct it!
      getLog() << Log::WARN << "This analysis uses the wrong jet definition: the " 
               << "paper just says 'a cone algorithm was applied to the CAL cells and jets "
               << "were reconstructed using the energies and positions of these cells'" << endl;
      addProjection(FastJets(fs, FastJets::KT, 0.7), "Jets");
    }


  // Book histograms
  void ZEUS_2001_S4815815::init() {
    /// @todo This doesn't seem to correspond to the plots in the paper (SPIRES 4730372)
    _histJetEt1 = bookHistogram1D("JetET1", "Jet transverse energy", 11, 14.0, 75.0);
  }


  // Do the analysis
  void ZEUS_2001_S4815815::analyze(const Event& event) {
    const FastJets& jets = applyProjection<FastJets>(event, "Jets");
    const size_t nj = jets.size();
    getLog() << Log::INFO << "Jet multiplicity = " << nj << endl;

    // Fill histograms
    PseudoJets jetList = jets.pseudoJets();
    for (PseudoJets::const_iterator j = jetList.begin(); j != jetList.end(); ++j) {
      _histJetEt1->fill(j->perp(), event.weight() );
    }
  }


  // Finalize
  void ZEUS_2001_S4815815::finalize() { }

}
