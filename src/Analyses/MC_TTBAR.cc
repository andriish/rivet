#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {
  
  class MC_TTBAR : public Analysis {

  public:

    MC_TTBAR()
      : Analysis("MC_TTBAR")
    {   }
    
    /// @name Analysis methods
    //@{
    void init() {
      addProjection(ChargedFinalState(-3.5, 3.5, 0.5*GeV), "CFS");
      addProjection(ChargedLeptons(ChargedFinalState(-3.5, 3.5, 30*GeV)), "LFS");
      addProjection(FastJets(FinalState(-2.5, 2.5, 0*GeV), FastJets::KT, 0.5), "JETS");

      _h_jet_0_pT = bookHistogram1D("jet_0_pT", 50, 0, 250);
      _h_jet_1_pT = bookHistogram1D("jet_1_pT", 50, 0, 250);
      _h_jet_2_pT = bookHistogram1D("jet_2_pT", 50, 0, 250);
      _h_jet_3_pT = bookHistogram1D("jet_3_pT", 50, 0, 250);

      _h_bjet_0_pT = bookHistogram1D("bjet_0_pT", 50, 0, 250);
      _h_bjet_1_pT = bookHistogram1D("bjet_1_pT", 50, 0, 250);

      _h_ljet_0_pT = bookHistogram1D("ljet_0_pT", 50, 0, 250);
      _h_ljet_1_pT = bookHistogram1D("ljet_1_pT", 50, 0, 250);

      _h_W_mass = bookHistogram1D("W_mass", 75, 30, 180);
      _h_t_mass = bookHistogram1D("t_mass", 150, 130, 430);
      _h_t_mass_W_cut = bookHistogram1D("t_mass_W_cut", 150, 130, 430);
      _h_W_comb_mass = bookHistogram1D("W_comb_mass", 75, 30, 180);
      _h_t_comb_mass = bookHistogram1D("t_comb_mass", 150, 130, 430);
    }
    
    void analyze(const Event& event) {
      double weight = event.weight();

      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      getLog() << Log::DEBUG << "Total charged multiplicity    = " 
               << cfs.size()  << endl;

      const ChargedLeptons& lfs = applyProjection<ChargedLeptons>(event, "LFS");
      getLog() << Log::DEBUG << "Charged lepton multiplicity   = " 
               << lfs.chargedLeptons().size()  << endl;
      if (lfs.chargedLeptons().size() != 1) {
        MSG_DEBUG("Event failed lepton cut");
        vetoEvent;
      }
      foreach (Particle lepton, lfs.chargedLeptons()) {
        getLog() << Log::DEBUG << "lepton pT = " << lepton.momentum().pT() << endl;
      }

      const FastJets& jetpro = applyProjection<FastJets>(event, "JETS");
      const Jets jets = jetpro.jetsByPt();
      getLog() << Log::DEBUG << "jet multiplicity = " << jets.size() << endl;

      if (jets.size() < 4) {
        getLog() << Log::DEBUG << "Event failed jet cut" << endl;
        vetoEvent;
      }

      _h_jet_0_pT->fill(jets[0].momentum().pT(), weight);
      _h_jet_1_pT->fill(jets[1].momentum().pT(), weight);
      _h_jet_2_pT->fill(jets[2].momentum().pT(), weight);
      _h_jet_3_pT->fill(jets[3].momentum().pT(), weight);

      if (jets[3].momentum().pT() < 35) {
        getLog() << Log::DEBUG << "Event failed jet cut" << endl;
        vetoEvent;
      }

      foreach (Jet jet, jets) {
        getLog() << Log::DEBUG << "jet pT = " << jet.momentum().pT() << endl;
      }

      Jets bjets, ljets;
      foreach (Jet jet, jets) {
        if (jet.momentum().pT() < 35*GeV) continue;
        if (jet.containsBottom())
          bjets.push_back(jet);
        else
          ljets.push_back(jet);
      }

      if (bjets.size() !=2) {
        getLog() << Log::DEBUG << "Event failed b-tagging cut" << endl;
        vetoEvent;
      }

      _h_bjet_0_pT->fill(bjets[0].momentum().pT(), weight);
      _h_bjet_1_pT->fill(bjets[1].momentum().pT(), weight);

      _h_ljet_0_pT->fill(ljets[0].momentum().pT(), weight);
      _h_ljet_1_pT->fill(ljets[1].momentum().pT(), weight);

      FourMomentum W  = ljets[0].momentum() + ljets[1].momentum();
      FourMomentum t1 = W + bjets[0].momentum();
      FourMomentum t2 = W + bjets[1].momentum();

      _h_W_mass->fill(mass(W), weight);
      _h_t_mass->fill(mass(t1), weight);
      _h_t_mass->fill(mass(t2), weight);
      if (mass(W) > 70 && mass(W) < 90) {
        getLog() << Log::INFO << "W found with mass " << W.mass() << endl;
        _h_t_mass_W_cut->fill(mass(t1), weight);
        _h_t_mass_W_cut->fill(mass(t2), weight);
      }

      _h_W_comb_mass->fill(mass(jets[0].momentum() + jets[1].momentum()), weight);
      _h_W_comb_mass->fill(mass(jets[0].momentum() + jets[2].momentum()), weight);
      _h_W_comb_mass->fill(mass(jets[0].momentum() + jets[3].momentum()), weight);
      _h_W_comb_mass->fill(mass(jets[1].momentum() + jets[2].momentum()), weight);
      _h_W_comb_mass->fill(mass(jets[1].momentum() + jets[3].momentum()), weight);
      _h_W_comb_mass->fill(mass(jets[2].momentum() + jets[3].momentum()), weight);

      _h_t_comb_mass->fill(mass(jets[0].momentum() + jets[1].momentum() + jets[2].momentum()), weight);
      _h_t_comb_mass->fill(mass(jets[0].momentum() + jets[1].momentum() + jets[3].momentum()), weight);
      _h_t_comb_mass->fill(mass(jets[0].momentum() + jets[2].momentum() + jets[3].momentum()), weight);
      _h_t_comb_mass->fill(mass(jets[1].momentum() + jets[2].momentum() + jets[3].momentum()), weight);
    }
    
    void finalize() {
      // No histos, so nothing to do!
    }
    //@}

  private:
    AIDA::IHistogram1D * _h_jet_0_pT;
    AIDA::IHistogram1D * _h_jet_1_pT;
    AIDA::IHistogram1D * _h_jet_2_pT;
    AIDA::IHistogram1D * _h_jet_3_pT;

    AIDA::IHistogram1D * _h_bjet_0_pT;
    AIDA::IHistogram1D * _h_bjet_1_pT;

    AIDA::IHistogram1D * _h_ljet_0_pT;
    AIDA::IHistogram1D * _h_ljet_1_pT;

    AIDA::IHistogram1D * _h_W_mass;
    AIDA::IHistogram1D * _h_t_mass;
    AIDA::IHistogram1D * _h_W_comb_mass;
    AIDA::IHistogram1D * _h_t_comb_mass;
    AIDA::IHistogram1D * _h_t_mass_W_cut;
  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_TTBAR> plugin_MC_TTBAR;

}
