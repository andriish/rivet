// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {

  /// @brief inclusive jet cross-section in pp at 200 GeV
  class STAR_2006_S6870392 : public Analysis {
  public:
    
    /// Constructor
    STAR_2006_S6870392()
      : Analysis("STAR_2006_S6870392")
    {
      setBeams(PROTON, PROTON);
      FinalState fs(-2.0, 2.0);
      addProjection(fs, "FS");
      // R=0.4, pTmin=0, seed_threshold=0.5:
      addProjection(FastJets(fs, FastJets::CDFMIDPOINT, 0.4, 0.0, 0.5), "MidpointJets");
    } 


    /// @name Analysis methods
    //@{ 

    /// Book histograms
    void init() {
      _h_jet_pT_MB = bookHistogram1D(1, 1, 1);
      _h_jet_pT_HT = bookHistogram1D(2, 1, 1);
    }

    /// Do the analysis 
    void analyze(const Event& event) {
      const double weight = event.weight();
      
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
    
    
    
    /// Finalize
    void finalize() {
      /// @todo Use the generator cross-section
      //_h_total_cross_section->fill(crossSection());
      normalize(_h_jet_pT_MB, 16603100);
      normalize(_h_jet_pT_HT, 1808234);
    }
    
    //@}
    
    
  private:
    
    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_jet_pT_MB;
    AIDA::IHistogram1D * _h_jet_pT_HT;
    //@}
    
  };
  
  
  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<STAR_2006_S6870392> plugin_STAR_2006_S6870392;
  
}
