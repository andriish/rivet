// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"

namespace Rivet {


  class EXAMPLE_SMEAR : public Analysis {
  public:

    /// Constructor
    EXAMPLE_SMEAR()
      : Analysis("EXAMPLE_SMEAR")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FastJets fj(FinalState(Cuts::abseta < 5), FastJets::ANTIKT, 0.4);
      addProjection(fj, "Jets0");

      SmearedJets sj1(fj, JET_SMEAR_IDENTITY);
      addProjection(sj1, "Jets1");

      SmearedJets sj2(fj, JET_SMEAR_ATLAS_RUN1,
                      [](const Jet& j){ return j.bTagged() ? 0.7*(1 - exp(-j.pT()/(10*GeV))) : 0.01; } );
      addProjection(sj2, "Jets2");

      SmearedJets sj3(fj,
                      [](const Jet& j){ return j; },
                      [](const Jet& j){ return j.bTagged() ? 0.7*(1 - exp(-j.pT()/(10*GeV))) : 0.01; },
                      JET_CTAG_PERFECT,
                      [](const Jet& j){ return 0.8; });
      addProjection(sj3, "Jets3");

      IdentifiedFinalState photons(Cuts::abseta < 5, PID::PHOTON);

      IdentifiedFinalState truthelectrons(Cuts::abseta < 5 && Cuts::pT > 10*GeV, {{PID::ELECTRON, PID::POSITRON}});
      addProjection(truthelectrons, "Electrons0");
      DressedLeptons dressedelectrons(photons, truthelectrons, 0.2);
      addProjection(dressedelectrons, "Electrons1");
      SmearedParticles recoelectrons(truthelectrons, ELECTRON_EFF_ATLAS_RUN1, ELECTRON_SMEAR_ATLAS_RUN1); //< @note Can't use dressedelectrons yet...
      addProjection(recoelectrons, "Electrons2");

      IdentifiedFinalState truthmuons(Cuts::abseta < 5 && Cuts::pT > 10*GeV, {{PID::MUON, PID::ANTIMUON}});
      addProjection(truthmuons, "Muons0");
      DressedLeptons dressedmuons(photons, truthmuons, 0.2);
      addProjection(dressedmuons, "Muons1");
      SmearedParticles recomuons(truthmuons, MUON_EFF_ATLAS_RUN1, MUON_SMEAR_ATLAS_RUN1); //< @note Can't use dressedmuons yet...
      addProjection(recomuons, "Muons2");

      TauFinder truthtaus(TauFinder::ANY, Cuts::abseta < 5 && Cuts::pT > 10*GeV);
      addProjection(truthtaus, "Taus0");
      DressedLeptons dressedtaus(photons, truthtaus, 0.2);
      addProjection(dressedtaus, "Taus1");
      SmearedParticles recotaus(truthtaus, TAU_EFF_ATLAS_RUN1, TAU_SMEAR_ATLAS_RUN1); //< @note Can't use dressedtaus yet...
      addProjection(recotaus, "Taus2");

      _h_nj_true = bookHisto1D("jet_N_true", 10, -0.5, 9.5);
      _h_nj_reco = bookHisto1D("jet_N_reco", 10, -0.5, 9.5);
      _h_j1pt_true = bookHisto1D("jet_pt1_true", 30, 0.0, 120);
      _h_j1pt_reco = bookHisto1D("jet_pt1_reco", 30, 0.0, 120);
      _h_j1eta_true = bookHisto1D("jet_eta1_true", 20, -5.0, 5.0);
      _h_j1eta_reco = bookHisto1D("jet_eta1_reco", 20, -5.0, 5.0);

      _h_ne_true = bookHisto1D("elec_N_true", 5, -0.5, 4.5);
      _h_ne_reco = bookHisto1D("elec_N_reco", 5, -0.5, 4.5);
      _h_e1pt_true = bookHisto1D("elec_pt1_true", 30, 0, 120);
      _h_e1pt_reco = bookHisto1D("elec_pt1_reco", 30, 0, 120);
      _h_e1eta_true = bookHisto1D("elec_eta1_true", 20, -5.0, 5.0);
      _h_e1eta_reco = bookHisto1D("elec_eta1_reco", 20, -5.0, 5.0);

      _h_nm_true = bookHisto1D("muon_N_true", 5, -0.5, 4.5);
      _h_nm_reco = bookHisto1D("muon_N_reco", 5, -0.5, 4.5);
      _h_m1pt_true = bookHisto1D("muon_pt1_true", 30, 0, 120);
      _h_m1pt_reco = bookHisto1D("muon_pt1_reco", 30, 0, 120);
      _h_m1eta_true = bookHisto1D("muon_eta1_true", 20, -5.0, 5.0);
      _h_m1eta_reco = bookHisto1D("muon_eta1_reco", 20, -5.0, 5.0);

      _h_nt_true = bookHisto1D("tau_N_true", 5, -0.5, 4.5);
      _h_nt_reco = bookHisto1D("tau_N_reco", 5, -0.5, 4.5);
      _h_t1pt_true = bookHisto1D("tau_pt1_true", 30, 0, 120);
      _h_t1pt_reco = bookHisto1D("tau_pt1_reco", 30, 0, 120);
      _h_t1eta_true = bookHisto1D("tau_eta1_true", 20, -5.0, 5.0);
      _h_t1eta_reco = bookHisto1D("tau_eta1_reco", 20, -5.0, 5.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const Jets jets0 = applyProjection<JetAlg>(event, "Jets0").jetsByPt(Cuts::pT > 10*GeV);
      const Jets jets1 = applyProjection<JetAlg>(event, "Jets1").jetsByPt(Cuts::pT > 10*GeV);
      const Jets jets2 = applyProjection<JetAlg>(event, "Jets2").jetsByPt(Cuts::pT > 10*GeV);
      const Jets jets3 = applyProjection<JetAlg>(event, "Jets3").jetsByPt(Cuts::pT > 10*GeV);
      MSG_DEBUG("Numbers of jets = " << jets0.size() << " true; "
               << jets1.size() << ", " << jets2.size() << ", " << jets3.size());
      _h_nj_true->fill(jets0.size(), weight);
      _h_nj_reco->fill(jets2.size(), weight);
      if (!jets0.empty()) {
        _h_j1pt_true->fill(jets0.front().pT()/GeV, weight);
        _h_j1eta_true->fill(jets0.front().eta(), weight);
      }
      if (!jets2.empty()) {
        _h_j1pt_reco->fill(jets2.front().pT()/GeV, weight);
        _h_j1eta_reco->fill(jets2.front().eta(), weight);
      }

      const Particles& elecs1 = applyProjection<ParticleFinder>(event, "Electrons1").particlesByPt();
      const Particles& elecs2 = applyProjection<ParticleFinder>(event, "Electrons2").particlesByPt();
      MSG_DEBUG("Numbers of electrons = " << elecs1.size() << " true; " << elecs2.size() << " reco");
      _h_ne_true->fill(elecs1.size(), weight);
      _h_ne_reco->fill(elecs2.size(), weight);
      if (!elecs1.empty()) {
        _h_e1pt_true->fill(elecs1.front().pT()/GeV, weight);
        _h_e1eta_true->fill(elecs1.front().eta(), weight);
      }
      if (!elecs2.empty()) {
        _h_e1pt_reco->fill(elecs2.front().pT()/GeV, weight);
        _h_e1eta_reco->fill(elecs2.front().eta(), weight);
      }

      const Particles& muons1 = applyProjection<ParticleFinder>(event, "Muons1").particlesByPt();
      const Particles& muons2 = applyProjection<ParticleFinder>(event, "Muons2").particlesByPt();
      MSG_DEBUG("Numbers of muons = " << muons1.size() << " true; " << muons2.size() << " reco");
      _h_nm_true->fill(muons1.size(), weight);
      _h_nm_reco->fill(muons2.size(), weight);
      if (!muons1.empty()) {
        _h_m1pt_true->fill(muons1.front().pT()/GeV, weight);
        _h_m1eta_true->fill(muons1.front().eta(), weight);
      }
      if (!muons2.empty()) {
        _h_m1pt_reco->fill(muons2.front().pT()/GeV, weight);
        _h_m1eta_reco->fill(muons2.front().eta(), weight);
      }

      const Particles& taus1 = applyProjection<ParticleFinder>(event, "Taus1").particlesByPt();
      const Particles& taus2 = applyProjection<ParticleFinder>(event, "Taus2").particlesByPt();
      MSG_DEBUG("Numbers of taus = " << taus1.size() << " true; " << taus2.size() << " reco");
      _h_nt_true->fill(taus1.size(), weight);
      _h_nt_reco->fill(taus2.size(), weight);
      if (!taus1.empty()) {
        _h_t1pt_true->fill(taus1.front().pT()/GeV, weight);
        _h_t1eta_true->fill(taus1.front().eta(), weight);
      }
      if (!taus2.empty()) {
        _h_t1pt_reco->fill(taus2.front().pT()/GeV, weight);
        _h_t1eta_reco->fill(taus2.front().eta(), weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_nj_true);
      normalize(_h_nj_reco);
      normalize(_h_j1pt_true, 1-_h_nj_true->bin(0).area());
      normalize(_h_j1pt_reco, 1-_h_nj_reco->bin(0).area());
      normalize(_h_j1eta_true, 1-_h_nj_true->bin(0).area());
      normalize(_h_j1eta_reco, 1-_h_nj_reco->bin(0).area());

      normalize(_h_ne_true);
      normalize(_h_ne_reco);
      normalize(_h_e1pt_true, 1-_h_ne_true->bin(0).area());
      normalize(_h_e1pt_reco, 1-_h_ne_reco->bin(0).area());
      normalize(_h_e1eta_true, 1-_h_ne_true->bin(0).area());
      normalize(_h_e1eta_reco, 1-_h_ne_reco->bin(0).area());

      normalize(_h_nm_true);
      normalize(_h_nm_reco);
      normalize(_h_m1pt_true, 1-_h_nm_true->bin(0).area());
      normalize(_h_m1pt_reco, 1-_h_nm_reco->bin(0).area());
      normalize(_h_m1eta_true, 1-_h_nm_true->bin(0).area());
      normalize(_h_m1eta_reco, 1-_h_nm_reco->bin(0).area());

      normalize(_h_nt_true);
      normalize(_h_nt_reco);
      normalize(_h_t1pt_true, 1-_h_nt_true->bin(0).area());
      normalize(_h_t1pt_reco, 1-_h_nt_reco->bin(0).area());
      normalize(_h_t1eta_true, 1-_h_nt_true->bin(0).area());
      normalize(_h_t1eta_reco, 1-_h_nt_reco->bin(0).area());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_nj_true, _h_nj_reco, _h_ne_true, _h_ne_reco,  _h_nm_true, _h_nm_reco,  _h_nt_true, _h_nt_reco;
    Histo1DPtr _h_j1pt_true, _h_j1pt_reco, _h_e1pt_true, _h_e1pt_reco,  _h_m1pt_true, _h_m1pt_reco,  _h_t1pt_true, _h_t1pt_reco;
    Histo1DPtr _h_j1eta_true, _h_j1eta_reco, _h_e1eta_true, _h_e1eta_reco,  _h_m1eta_true, _h_m1eta_reco,  _h_t1eta_true, _h_t1eta_reco;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EXAMPLE_SMEAR);


}
