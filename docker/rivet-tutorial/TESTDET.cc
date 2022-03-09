// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/Smearing.hh"

namespace Rivet {


  class TESTDET : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TESTDET);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState calofs(Cuts::abseta < 4.9);

      FastJets fj(calofs, FastJets::ANTIKT, 0.4);
      declare(fj, "TrueJets");
      SmearedJets sj(fj, JET_SMEAR_ATLAS_RUN2);
      declare(sj, "RecoJets");

      MissingMomentum mm(calofs);
      declare(mm, "TrueMET");
      SmearedMET sm(mm, MET_SMEAR_ATLAS_RUN2);
      declare(sm, "RecoMET");

      FinalState truthelectrons(Cuts::abseta < 2.5 && Cuts::pT > 10*GeV && Cuts::abspid == PID::ELECTRON);
      declare(truthelectrons, "TrueElecs");
      SmearedParticles recoelectrons(truthelectrons, ELECTRON_EFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2);
      declare(recoelectrons, "RecoElecs");

      FinalState truthmuons(Cuts::abseta < 2.7 && Cuts::pT > 10*GeV && Cuts::abspid == PID::ELECTRON);
      declare(truthmuons, "TrueMuons");
      SmearedParticles recomuons(truthmuons, MUON_EFF_ATLAS_RUN1, MUON_SMEAR_ATLAS_RUN1);
      declare(recomuons, "RecoMuons");


      _h_nj_true = bookHisto1D("jet_N_true", 10, -0.5, 9.5);
      _h_nj_reco = bookHisto1D("jet_N_reco", 10, -0.5, 9.5);
      _h_nbj_true = bookHisto1D("jet_Nb_true", 10, -0.5, 9.5);
      _h_nbj_reco = bookHisto1D("jet_Nb_reco", 10, -0.5, 9.5);
      _h_j1pt_true = bookHisto1D("jet_pt1_true", logspace(30, 1, 400));
      _h_j1pt_reco = bookHisto1D("jet_pt1_reco", logspace(30, 1, 400));
      _h_j1eta_true = bookHisto1D("jet_eta1_true", 20, -5.0, 5.0);
      _h_j1eta_reco = bookHisto1D("jet_eta1_reco", 20, -5.0, 5.0);

      _h_ne_true = bookHisto1D("elec_N_true", 5, -0.5, 4.5);
      _h_ne_reco = bookHisto1D("elec_N_reco", 5, -0.5, 4.5);
      _h_e1pt_true = bookHisto1D("elec_pt1_true", logspace(30, 1, 400));
      _h_e1pt_reco = bookHisto1D("elec_pt1_reco", logspace(30, 1, 400));
      _h_e1eta_true = bookHisto1D("elec_eta1_true", 20, -5.0, 5.0);
      _h_e1eta_reco = bookHisto1D("elec_eta1_reco", 20, -5.0, 5.0);

      _h_nm_true = bookHisto1D("muon_N_true", 5, -0.5, 4.5);
      _h_nm_reco = bookHisto1D("muon_N_reco", 5, -0.5, 4.5);
      _h_m1pt_true = bookHisto1D("muon_pt1_true", logspace(30, 1, 400));
      _h_m1pt_reco = bookHisto1D("muon_pt1_reco", logspace(30, 1, 400));
      _h_m1eta_true = bookHisto1D("muon_eta1_true", 20, -5.0, 5.0);
      _h_m1eta_reco = bookHisto1D("muon_eta1_reco", 20, -5.0, 5.0);

      _h_met_true = bookHisto1D("met_true", logspace(30, 1, 400));
      _h_met_reco = bookHisto1D("met_reco", logspace(30, 1, 400));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const Jets tjets = apply<JetAlg>(event, "TrueJets").jetsByPt(Cuts::pT > 40*GeV);
      const Jets rjets = apply<JetAlg>(event, "RecoJets").jetsByPt(Cuts::pT > 40*GeV);
      MSG_DEBUG("Numbers of jets = " << tjets.size() << " true; " << rjets.size() << " reco");
      _h_nj_true->fill(tjets.size(), weight);
      _h_nj_reco->fill(rjets.size(), weight);
      const Jets tbjets = filter_select(tjets, [](const Jet& j){ return j.bTagged(Cuts::pT > 5*GeV); });
      const Jets rbjets = filter_select(rjets, [](const Jet& j){ return j.bTagged(Cuts::pT > 5*GeV); });
      MSG_DEBUG("Numbers of b-jets = " << tbjets.size() << " true; " << rbjets.size() << " reco");
      _h_nbj_true->fill(tbjets.size(), weight);
      _h_nbj_reco->fill(rbjets.size(), weight);
      if (!tjets.empty()) {
        _h_j1pt_true->fill(tjets.front().pT()/GeV, weight);
        _h_j1eta_true->fill(tjets.front().eta(), weight);
      }
      if (!rjets.empty()) {
        _h_j1pt_reco->fill(rjets.front().pT()/GeV, weight);
        _h_j1eta_reco->fill(rjets.front().eta(), weight);
      }

      const Particles& telecs = apply<ParticleFinder>(event, "TrueElecs").particlesByPt(20*GeV);
      Particles relecs = apply<ParticleFinder>(event, "RecoElecs").particlesByPt(20*GeV);
      ifilter_select(relecs, particleEffFilter(ELECTRON_IDEFF_ATLAS_RUN2_TIGHT));
      MSG_DEBUG("Numbers of electrons = " << telecs.size() << " true; " << relecs.size() << " reco");
      _h_ne_true->fill(telecs.size(), weight);
      _h_ne_reco->fill(relecs.size(), weight);
      if (!telecs.empty()) {
        _h_e1pt_true->fill(telecs.front().pT()/GeV, weight);
        _h_e1eta_true->fill(telecs.front().eta(), weight);
      }
      if (!relecs.empty()) {
        _h_e1pt_reco->fill(relecs.front().pT()/GeV, weight);
        _h_e1eta_reco->fill(relecs.front().eta(), weight);
      }

      const Particles& tmuons = apply<ParticleFinder>(event, "TrueMuons").particlesByPt(20*GeV);
      const Particles& rmuons = apply<ParticleFinder>(event, "RecoMuons").particlesByPt(20*GeV);
      MSG_DEBUG("Numbers of muons = " << tmuons.size() << " true; " << rmuons.size() << " reco");
      _h_nm_true->fill(tmuons.size(), weight);
      _h_nm_reco->fill(rmuons.size(), weight);
      if (!tmuons.empty()) {
        _h_m1pt_true->fill(tmuons.front().pT()/GeV, weight);
        _h_m1eta_true->fill(tmuons.front().eta(), weight);
      }
      if (!rmuons.empty()) {
        _h_m1pt_reco->fill(rmuons.front().pT()/GeV, weight);
        _h_m1eta_reco->fill(rmuons.front().eta(), weight);
      }

      const double tmet = apply<MissingMomentum>(event, "TrueMET").met();
      const double rmet = apply<SmearedMET>(event, "RecoMET").met();
      _h_met_true->fill(tmet/GeV, weight);
      _h_met_reco->fill(rmet/GeV, weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize(_h_nj_true);
      // normalize(_h_nj_reco);
      // normalize(_h_nbj_true);
      // normalize(_h_nbj_reco);
      // normalize(_h_j1pt_true, 1-_h_nj_true->bin(0).area());
      // normalize(_h_j1pt_reco, 1-_h_nj_reco->bin(0).area());
      // normalize(_h_j1eta_true, 1-_h_nj_true->bin(0).area());
      // normalize(_h_j1eta_reco, 1-_h_nj_reco->bin(0).area());

      // normalize(_h_ne_true);
      // normalize(_h_ne_reco);
      // normalize(_h_e1pt_true, 1-_h_ne_true->bin(0).area());
      // normalize(_h_e1pt_reco, 1-_h_ne_reco->bin(0).area());
      // normalize(_h_e1eta_true, 1-_h_ne_true->bin(0).area());
      // normalize(_h_e1eta_reco, 1-_h_ne_reco->bin(0).area());

      // normalize(_h_nm_true);
      // normalize(_h_nm_reco);
      // normalize(_h_m1pt_true, 1-_h_nm_true->bin(0).area());
      // normalize(_h_m1pt_reco, 1-_h_nm_reco->bin(0).area());
      // normalize(_h_m1eta_true, 1-_h_nm_true->bin(0).area());
      // normalize(_h_m1eta_reco, 1-_h_nm_reco->bin(0).area());

      // normalize(_h_met_true);
      // normalize(_h_met_reco);

    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_nj_true, _h_nj_reco, _h_nbj_true, _h_nbj_reco, _h_ne_true, _h_ne_reco,  _h_nm_true, _h_nm_reco;
    Histo1DPtr _h_j1pt_true, _h_j1pt_reco, _h_e1pt_true, _h_e1pt_reco,  _h_m1pt_true, _h_m1pt_reco;
    Histo1DPtr _h_j1eta_true, _h_j1eta_reco, _h_e1eta_true, _h_e1eta_reco,  _h_m1eta_true, _h_m1eta_reco;
    Histo1DPtr _h_met_true, _h_met_reco;
    /// @}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TESTDET);


}
