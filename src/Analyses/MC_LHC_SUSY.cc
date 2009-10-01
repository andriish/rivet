// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  /* Basic SUSY type validation analysis for the LHC
   * @author Andy Buckley
   */ 
  class MC_LHC_SUSY : public Analysis {
  public:
    
    /// Constructor
    MC_LHC_SUSY()
      : Analysis("MC_LHC_SUSY")
    { 
      setBeams(PROTON, PROTON);
    }
    
    
    /// @name Analysis methods
    //@{
    
    // Book histograms
    void init() {
      // Basic final state
      const FinalState fs(-4.0, 4.0, 0.5*GeV);

      // Tracks and jets
      addProjection(ChargedFinalState(fs), "Tracks");
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.7), "Jets");

      IdentifiedFinalState photonfs(fs);
      photonfs.acceptId(PHOTON);
      addProjection(photonfs, "AllPhotons");
      /// @todo Isolated photons

      IdentifiedFinalState efs(fs);
      efs.acceptIdPair(ELECTRON);
      addProjection(efs, "Electrons");

      IdentifiedFinalState mufs(fs);
      efs.acceptIdPair(MUON);
      addProjection(mufs, "Muons");

      addProjection(TotalVisibleMomentum(fs), "MET");
      LeadingParticlesFinalState lpfs(fs);
      lpfs.addParticleIdPair(ELECTRON);
      lpfs.addParticleIdPair(MUON);
      addProjection(lpfs, "LeadingParticles");

      _hist_n_trk   = bookHistogram1D("n-trk", 50, 0.5, 300.5);
      _hist_phi_trk = bookHistogram1D("phi-trk", 50, -PI, PI);
      _hist_eta_trk = bookHistogram1D("eta-trk", 50, -4, 4);
      _hist_pt_trk  = bookHistogram1D("pt-trk", 100, 0.0, 1500);

      _hist_n_jet   = bookHistogram1D("n-jet", 21, -0.5, 20.5);
      _hist_phi_jet = bookHistogram1D("phi-jet", 50, -PI, PI);
      _hist_eta_jet = bookHistogram1D("eta-jet", 50, -4, 4);
      _hist_pt_jet  = bookHistogram1D("pt-jet", 100, 0.0, 1500);

      _hist_n_e   = bookHistogram1D("n-e", 11, -0.5, 10.5);
      _hist_phi_e = bookHistogram1D("phi-e", 50, -PI, PI);
      _hist_eta_e = bookHistogram1D("eta-e", 50, -4, 4);
      _hist_pt_e  = bookHistogram1D("pt-e", 100, 0.0, 500);

      _hist_n_mu   = bookHistogram1D("n-mu", 11, -0.5, 10.5);
      _hist_phi_mu = bookHistogram1D("phi-mu", 50, -PI, PI);
      _hist_eta_mu = bookHistogram1D("eta-mu", 50, -4, 4);
      _hist_pt_mu  = bookHistogram1D("pt-mu", 100, 0.0, 500);

      _hist_n_gamma   = bookHistogram1D("n-gamma", 11, -0.5, 10.5);
      _hist_phi_gamma = bookHistogram1D("phi-gamma", 50, -PI, PI);
      _hist_eta_gamma = bookHistogram1D("eta-gamma", 50, -4, 4);
      _hist_pt_gamma  = bookHistogram1D("pt-gamma", 100, 0.0, 500);

      /// @todo Isolated photons

      _hist_met = bookHistogram1D("Etmiss", 100, 0.0, 1500);

      _hist_mll_ossf_ee   = bookHistogram1D("mll-ossf-ee", 50, 0.0, 500);
      _hist_mll_ossf_mumu = bookHistogram1D("mll-ossf-mumu", 50, 0.0, 500);
      _hist_mll_osof_emu  = bookHistogram1D("mll-osof-mumu", 50, 0.0, 500);

      // LSP eta, pT, phi, mass
    }


    // Do the analysis
    void analyze(const Event& evt) {
      const FinalState& tracks = applyProjection<FinalState>(evt, "Tracks");
      if (tracks.particles().empty()) {
        getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
        vetoEvent;
      }

      // Get event weight
      const double weight = evt.weight();

      // Fill track histos
      _hist_n_trk->fill(tracks.size(), weight);
      foreach (const Particle& t, tracks.particles()) {
        const FourMomentum& p = t.momentum();
        _hist_phi_trk->fill(p.phi(), weight);
        _hist_eta_trk->fill(p.eta(), weight);
        _hist_pt_trk->fill(p.pT()/GeV, weight);
      }

      // Get jets and fill jet histos
      const FastJets& jetpro = applyProjection<FastJets>(evt, "Jets");
      const Jets jets = jetpro.jetsByPt();
      getLog() << Log::DEBUG << "Jet multiplicity = " << jets.size() << endl;
      _hist_n_jet->fill(jets.size(), weight);
      foreach (const Jet& j, jets) {
        const FourMomentum& pj = j.momentum();
        _hist_phi_jet->fill(pj.phi(), weight);
        _hist_eta_jet->fill(pj.eta(), weight);
        _hist_pt_jet->fill(pj.pT()/GeV, weight);
      }

      /// @todo Resum photons around electrons

      // Fill final state electron/positron histos
      const FinalState& efs = applyProjection<FinalState>(evt, "Electrons");
      _hist_n_e->fill(efs.size(), weight);
      foreach (const Particle& e, efs.particles()) {
        const FourMomentum& p = e.momentum();
        _hist_phi_e->fill(p.phi(), weight);
        _hist_eta_e->fill(p.eta(), weight);
        _hist_pt_e->fill(p.pT()/GeV, weight);
      }

      /// @todo Resum photons around muons

      // Fill final state muon/antimuon histos
      const FinalState& mufs = applyProjection<FinalState>(evt, "Muons");
      _hist_n_e->fill(mufs.size(), weight);
      foreach (const Particle& mu, efs.particles()) {
        const FourMomentum& p = mu.momentum();
        _hist_phi_mu->fill(p.phi(), weight);
        _hist_eta_mu->fill(p.eta(), weight);
        _hist_pt_mu->fill(p.pT()/GeV, weight);
      }

      // Fill final state non-isolated photon histos
      const FinalState& allphotonfs = applyProjection<FinalState>(evt, "AllPhotons");
      _hist_n_e->fill(allphotonfs.size(), weight);
      foreach (const Particle& ph, allphotonfs.particles()) {
        const FourMomentum& p = ph.momentum();
        _hist_phi_mu->fill(p.phi(), weight);
        _hist_eta_mu->fill(p.eta(), weight);
        _hist_pt_mu->fill(p.pT()/GeV, weight);
      }

      /// @todo Isolated photons

      // Calculate and fill missing Et histos
      const TotalVisibleMomentum& met = applyProjection<TotalVisibleMomentum>(evt, "MET");
      _hist_met->fill(met.scalarET()/GeV);

      // Choose highest-pT leptons of each sign and flavour for dilepton mass edges
      const FinalState& lpfs = applyProjection<FinalState>(evt, "LeadingParticles");
      bool eplus_ok, eminus_ok, muplus_ok, muminus_ok;
      FourMomentum peplus, peminus, pmuplus, pmuminus;
      foreach (const Particle& p, lpfs.particles()) {
        const PdgId pid = p.pdgId();
        if (pid == ELECTRON) {
          eminus_ok = true;
          peminus = p.momentum();
        } else if (pid == POSITRON) {
          eplus_ok = true;
          peplus = p.momentum();
        } else if (pid == MUON) {
          muminus_ok = true;
          pmuminus = p.momentum();
        } else if (pid == ANTIMUON) {
          muplus_ok = true;
          pmuplus = p.momentum();
        } else {
          throw Error("Unexpected particle type in leading particles FS!");
        }
      }
      // m_ee
      if (eminus_ok && eplus_ok) {
        const double m_ee = FourMomentum(peplus + peminus).mass();
        _hist_mll_ossf_ee->fill(m_ee/GeV, weight);
      }
      // m_mumu
      if (muminus_ok && muplus_ok) {
        const double m_mumu = FourMomentum(pmuplus + pmuminus).mass();
        _hist_mll_ossf_mumu->fill(m_mumu/GeV, weight);
      }
      // m_emu (both configurations)
      if (eminus_ok && muplus_ok) {
        const double m_emu = FourMomentum(pmuplus + peminus).mass();
        _hist_mll_osof_emu->fill(m_emu/GeV, weight);
      }
      if (muminus_ok && eplus_ok) {
        const double m_mue = FourMomentum(peplus + pmuminus).mass();
        _hist_mll_osof_emu->fill(m_mue/GeV, weight);
      }

    }
    
    
    void finalize() {  
      /// @todo Normalisations
    }

    //@}    
    

  private:
    
    AIDA::IHistogram1D *_hist_n_trk, *_hist_phi_trk, *_hist_eta_trk, *_hist_pt_trk;
    AIDA::IHistogram1D *_hist_n_jet, *_hist_phi_jet, *_hist_eta_jet, *_hist_pt_jet;
    AIDA::IHistogram1D *_hist_n_e, *_hist_phi_e, *_hist_eta_e, *_hist_pt_e;
    AIDA::IHistogram1D *_hist_n_mu, *_hist_phi_mu, *_hist_eta_mu, *_hist_pt_mu;
    AIDA::IHistogram1D *_hist_n_gamma, *_hist_phi_gamma, *_hist_eta_gamma, *_hist_pt_gamma;
    /// @todo Isolated photons
    AIDA::IHistogram1D *_hist_met;
    AIDA::IHistogram1D *_hist_mll_ossf_ee, *_hist_mll_ossf_mumu, *_hist_mll_osof_emu;
    
  };
  
  
  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_LHC_SUSY> plugin_MC_LHC_SUSY;
  
}
