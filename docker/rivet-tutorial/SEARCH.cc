// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/Smearing.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// Rivet search-analysis demo
  class SEARCH : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(SEARCH);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise projections
      FinalState fs(Cuts::abseta < 4.9);
      FastJets tjets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::DECAY, JetAlg::Invisibles::NONE);
      FinalState telecs(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.5);
      FinalState tmuons(Cuts::abspid == PID::MUON && Cuts::abseta < 2.5);
      MissingMomentum tmet(fs);
      //
      SmearedJets rjets(tjets, JET_SMEAR_CMS_RUN2, JET_BTAG_EFFS(0.7, 0.1, 1/120.));
      SmearedParticles relecs(telecs,
                              [](const Particle& e){
                                return 0.85 * (1-sqr(e.eta()/5)) * (1 - 0.1*exp(10 - e.pT()/2/GeV));
                              }, ELECTRON_SMEAR_CMS_RUN2);
      SmearedParticles rmuons(tmuons, PARTICLE_EFF_CONST(0.8), MUON_SMEAR_CMS_RUN2);
      SmearedMET rmet(tmet, MET_SMEAR_CMS_RUN2);

      // Register truth or reco projections according to the RECO analysis option
      bool reco = (getOption("RECO") != "0");
      if (reco) {
        MSG_INFO("Using detector smearing");
        declare(rjets, "Jets");
        declare(relecs, "Elecs");
        declare(rmuons, "Muons");
        declare(rmet, "MET");
      } else {
        declare(tjets, "Jets");
        declare(telecs, "Elecs");
        declare(tmuons, "Muons");
        declare(tmet, "MET");
      }

      // Book histograms
      book(_hists["j_pt"], "j_pt", logspace(20, 30, 500));
      book(_hists["j_eta"], "j_eta", linspace(20, -5, 5));
      book(_hists["jb_pt"], "jb_pt", logspace(20, 30, 500));
      book(_hists["jb_eta"], "jb_eta", linspace(20, -5, 5));
      book(_hists["le_pt"], "le_pt", logspace(20, 30, 500));
      book(_hists["le_eta"], "le_eta", linspace(20, -2.5, 2.5));
      book(_hists["lm_pt"], "lm_pt", logspace(20, 30, 500));
      book(_hists["lm_eta"], "lm_eta", linspace(20, -2.5, 2.5));

      // Book yield counters
      book(_srs["SR1L"], "SR1L");
      book(_srs["SR2L"], "SR2L");
      book(_srs["SR3L"], "SR3L");

      // Book cutflows
      strings cutnames = {"3 jets", "2 bjets", "1 lep", "HT", "MET", "SR"};
      _flows.addCutflow("SR1L", cutnames);
      _flows.addCutflow("SR2L", cutnames);
      _flows.addCutflow("SR3L", cutnames);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      _flows.fillinit();

      // Retrieve jets and leptons
      Jets jets = apply<JetFinder>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 4);
      Particles elecs = apply<ParticleFinder>(event, "Elecs").particlesByPt(Cuts::pT > 20*GeV);
      Particles muons = apply<ParticleFinder>(event, "Muons").particlesByPt(Cuts::pT > 20*GeV);

      // Remove all jets within dR < 0.2 of an electron, then electrons within 0.4 of a jet, then muons
      idiscardIfAnyDeltaRLess(jets, elecs, 0.2);
      idiscardIfAnyDeltaRLess(elecs, jets, 0.4);
      idiscardIfAnyDeltaRLess(muons, filter_select(jets, [](const Jet& j){
                                       return j.particles(Cuts::abscharge > 0).size() > 4;
                                     }), 0.4);

      // Get b-jets
      Jets bjets = filter_select(jets, hasBTag(Cuts::abseta < 2.5));

      // Fill histograms
      for (const Jet& j : jets) {
        _hists["j_pt"]->fill(j.pT()/GeV);
        _hists["j_eta"]->fill(j.eta());
      }
      for (const Jet& j : bjets) {
        _hists["jb_pt"]->fill(j.pT()/GeV);
        _hists["jb_eta"]->fill(j.eta());
      }
      for (const Particle& e : elecs) {
        _hists["le_pt"]->fill(e.pT()/GeV);
        _hists["le_eta"]->fill(e.eta());
      }
      for (const Particle& m : muons) {
        _hists["lm_pt"]->fill(m.pT()/GeV);
        _hists["lm_eta"]->fill(m.eta());
      }

      // Compute HT and MET
      double ht = sum(jets, Kin::pT, 0.0);
      double met = apply<METFinder>(event, "MET").missingPt();
      size_t nleps = elecs.size() + muons.size();

      // Event selection
      if (!_flows.fillnext(jets.size() >= 3)) vetoEvent;
      if (!_flows.fillnext(count(bjets, pTGtr(60*GeV)) >= 2)) vetoEvent;
      if (!_flows.fillnext(nleps > 0)) vetoEvent;
      if (!_flows.fillnext(ht > 800*GeV)) vetoEvent;
      if (!_flows.fillnext(met > 200*GeV)) vetoEvent;

      // Fill signal regions
      size_t lepindex = (nleps < 3) ? nleps-1 : 2;
      string srname = "SR" + to_string(lepindex+1) + "L";
      _flows[srname].fillnext(true);
      _srs[srname]->fill();

    }


    /// Normalise SR counters to give the expected event yields
    void finalize() {
      MSG_INFO(_flows);
      scale(_srs, 139*crossSection()/femtobarn/sumW());
    }

    /// @}

    /// Histograms
    map<string, Histo1DPtr> _hists;

    /// Signal region yield counters
    map<string, CounterPtr> _srs;

    /// Cut-flow monitors
    Cutflows _flows;

  };


  DECLARE_RIVET_PLUGIN(SEARCH);

}
