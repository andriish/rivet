#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {


  class MC_TTBAR : public Analysis {
  public:

    /// Minimal constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_TTBAR);


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {

      _mode = 1; string pre = "onelep_"; // default is single-lepton decay mode
      if ( getOption("TTMODE") == "ALLHAD" ) { _mode = 0; pre = "allhad_"; }
      if ( getOption("TTMODE") == "ONELEP" ) { _mode = 1; pre = "onelep_"; }
      if ( getOption("TTMODE") == "TWOLEP" ) { _mode = 2; pre = "twolep_"; }
      if ( getOption("TTMODE") == "ANYLEP" ) { _mode = 3; pre = "anylep_"; }

      // A FinalState is used to select particles within |eta| < 4.2 and with pT
      // > 30 GeV, out of which the ChargedLeptons projection picks only the
      // electrons and muons, to be accessed later as "LFS".
      ChargedLeptons lfs(FinalState(Cuts::abseta < 4.2 && Cuts::pT > 30*GeV));
      declare(lfs, "LFS");

      // A second FinalState is used to select all particles in |eta| < 4.2,
      // with no pT cut. This is used to construct jets and measure missing
      // transverse energy.
      VetoedFinalState fs(FinalState(Cuts::abseta < 4.2));
      fs.addVetoOnThisFinalState(lfs);
      declare(FastJets(fs, FastJets::ANTIKT, 0.6), "Jets");
      declare(MissingMomentum(fs), "MissingET");

      // Booking of histograms
      _h["njets"] = bookHisto1D(pre + "jet_mult", 11, -0.5, 10.5);
      //
      _h["jet_1_pT"] = bookHisto1D(pre + "jet_1_pT", logspace(50, 20.0, 500.0));
      _h["jet_2_pT"] = bookHisto1D(pre + "jet_2_pT", logspace(50, 20.0, 400.0));
      _h["jet_3_pT"] = bookHisto1D(pre + "jet_3_pT", logspace(50, 20.0, 300.0));
      _h["jet_4_pT"] = bookHisto1D(pre + "jet_4_pT", logspace(50, 20.0, 200.0));
      _h["jet_HT"]   = bookHisto1D(pre + "jet_HT", logspace(50, 100.0, 2000.0));
      //
      _h["bjet_1_pT"] = bookHisto1D(pre + "jetb_1_pT", logspace(50, 20.0, 400.0));
      _h["bjet_2_pT"] = bookHisto1D(pre + "jetb_2_pT", logspace(50, 20.0, 300.0));
      //
      _h["ljet_1_pT"] = bookHisto1D(pre + "jetl_1_pT", logspace(50, 20.0, 400.0));
      _h["ljet_2_pT"] = bookHisto1D(pre + "jetl_2_pT", logspace(50, 20.0, 300.0));
      //
      _h["W_mass"] = bookHisto1D(pre + "W_mass", 75, 30, 180);
      _h["t_mass"] = bookHisto1D(pre + "t_mass", 150, 130, 430);
      _h["t_mass_W_cut"] = bookHisto1D(pre + "t_mass_W_cut", 150, 130, 430);
      //
      _h["jetb_1_jetb_2_dR"]   = bookHisto1D(pre + "jetb_1_jetb_2_dR", 20, 0.0, 7.0);
      _h["jetb_1_jetb_2_deta"] = bookHisto1D(pre + "jetb_1_jetb_2_deta", 20, 0.0, 7.0);
      _h["jetb_1_jetb_2_dphi"] = bookHisto1D(pre + "jetb_1_jetb_2_dphi", 20, 0.0, M_PI);
      _h["jetb_1_jetl_1_dR"]   = bookHisto1D(pre + "jetb_1_jetl_1_dR", 20, 0.0, 7.0);
      _h["jetb_1_jetl_1_deta"] = bookHisto1D(pre + "jetb_1_jetl_1_deta", 20, 0.0, 7.0);
      _h["jetb_1_jetl_1_dphi"] = bookHisto1D(pre + "jetb_1_jetl_1_dphi", 20, 0.0, M_PI);
      _h["jetl_1_jetl_2_dR"]   = bookHisto1D(pre + "jetl_1_jetl_2_dR", 20, 0.0, 7.0);
      _h["jetl_1_jetl_2_deta"] = bookHisto1D(pre + "jetl_1_jetl_2_deta", 20, 0.0, 7.0);
      _h["jetl_1_jetl_2_dphi"] = bookHisto1D(pre + "jetl_1_jetl_2_dphi", 20, 0.0, M_PI);
      _h["jetb_1_W_dR"]        = bookHisto1D(pre + "jetb_1_W_dR", 20, 0.0, 7.0);
      _h["jetb_1_W_deta"]      = bookHisto1D(pre + "jetb_1_W_deta", 20, 0.0, 7.0);
      _h["jetb_1_W_dphi"]      = bookHisto1D(pre + "jetb_1_W_dphi", 20, 0.0, M_PI);
      if (_mode > 0) {
        _h["jetb_1_l_dR"]        = bookHisto1D(pre + "jetb_1_l_dR", 20, 0.0, 7.0);
        _h["jetb_1_l_deta"]      = bookHisto1D(pre + "jetb_1_l_deta", 20, 0.0, 7.0);
        _h["jetb_1_l_dphi"]      = bookHisto1D(pre + "jetb_1_l_dphi", 20, 0.0, M_PI);
        _h["jetb_1_l_mass"]      = bookHisto1D(pre + "jetb_1_l_mass", 40, 0.0, 500.0);
        if (_mode > 1) {
          _h["jetb_1_l2_dR"]       = bookHisto1D(pre + "jetb_1_l2_dR", 20, 0.0, 7.0);
          _h["jetb_1_l2_deta"]     = bookHisto1D(pre + "jetb_1_l2_deta", 20, 0.0, 7.0);
          _h["jetb_1_l2_dphi"]     = bookHisto1D(pre + "jetb_1_l2_dphi", 20, 0.0, M_PI);
          _h["jetb_1_l2_mass"]     = bookHisto1D(pre + "jetb_1_l2_mass", 40, 0.0, 500.0);
        }
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      // Use the "LFS" projection to require at least one hard charged
      // lepton. This is an experimental signature for the leptonically decaying
      // W. This helps to reduce pure QCD backgrounds.
      const ChargedLeptons& lfs = apply<ChargedLeptons>(event, "LFS");
      MSG_DEBUG("Charged lepton multiplicity = " << lfs.chargedLeptons().size());
      for (const Particle& lepton : lfs.chargedLeptons()) {
        MSG_DEBUG("Lepton pT = " << lepton.pT());
      }

      size_t nLeps = lfs.chargedLeptons().size();
      bool leptonMultiFail = _mode == 3 && nLeps == 0; // non-all-hadronic
      leptonMultiFail |= _mode == 2 && nLeps != 2; // dilepton
      leptonMultiFail |= _mode == 1 && nLeps != 1; // single lepton
      leptonMultiFail |= _mode == 0 && nLeps != 0; // all-hadronic
      if (leptonMultiFail) {
        MSG_DEBUG("Event failed lepton multiplicity cut");
        vetoEvent;
      }

      // Use a missing ET cut to bias toward events with a hard neutrino from
      // the leptonically decaying W. This helps to reduce pure QCD backgrounds.
      // not applied in all-hadronic mode
      const MissingMomentum& met = apply<MissingMomentum>(event, "MissingET");
      MSG_DEBUG("Vector ET = " << met.vectorEt().mod() << " GeV");
      if (_mode > 0 && met.vectorEt().mod() < 30*GeV) {
        MSG_DEBUG("Event failed missing ET cut");
        vetoEvent;
      }

      // Use the "Jets" projection to check that there are at least 4 jets of
      // any pT. Getting the jets sorted by pT ensures that the first jet is the
      // hardest, and so on. We apply no pT cut here only because we want to
      // plot all jet pTs to help optimise our jet pT cut.
      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets alljets = jetpro.jetsByPt();
      if (alljets.size() < 4) {
        MSG_DEBUG("Event failed jet multiplicity cut");
        vetoEvent;
      }

      // Update passed-cuts counter and fill all-jets histograms
      _h["jet_1_pT"]->fill(alljets[0].pT()/GeV, weight);
      _h["jet_2_pT"]->fill(alljets[1].pT()/GeV, weight);
      _h["jet_3_pT"]->fill(alljets[2].pT()/GeV, weight);
      _h["jet_4_pT"]->fill(alljets[3].pT()/GeV, weight);

      // Insist that the hardest 4 jets pass pT hardness cuts. If we don't find
      // at least 4 such jets, we abandon this event.
      const Jets jets = jetpro.jetsByPt(30*GeV);
      _h["njets"]->fill(jets.size(), weight);
      double ht = 0.0;
      for (const Jet& j : jets) { ht += j.pT(); }
      _h["jet_HT"]->fill(ht/GeV, weight);
      if (jets.size() < 4 ||
          jets[0].pT() < 60*GeV ||
          jets[1].pT() < 50*GeV ||
          jets[3].pT() < 30*GeV) {
        MSG_DEBUG("Event failed jet cuts");
        vetoEvent;
      }

      // Sort the jets into b-jets and light jets. We expect one hard b-jet from
      // each top decay, so our 4 hardest jets should include two b-jets. The
      // Jet::bTagged() method is equivalent to perfect experimental
      // b-tagging, in a generator-independent way.
      Jets bjets, ljets;
      for (const Jet& jet : jets) {
        // // Don't count jets that overlap with the hard leptons
        bool isolated = true;
        for (const Particle& lepton : lfs.chargedLeptons()) {
          if (deltaR(jet.momentum(), lepton.momentum()) < 0.3) {
            isolated = false;
            break;
          }
        }
        if (!isolated) {
          MSG_DEBUG("Jet failed lepton isolation cut");
          break;
        }
        if (jet.bTagged()) {
          bjets.push_back(jet);
        } else {
          ljets.push_back(jet);
        }
      }
      MSG_DEBUG("Number of b-jets = " << bjets.size());
      MSG_DEBUG("Number of l-jets = " << ljets.size());
      if (bjets.size() != 2) {
        MSG_DEBUG("Event failed post-lepton-isolation b-tagging cut");
        vetoEvent;
      }
      if (ljets.size() < 2) {
        MSG_DEBUG("Event failed since not enough light jets remaining after lepton-isolation");
        vetoEvent;
      }

      // Plot the pTs of the identified jets.
      _h["bjet_1_pT"]->fill(bjets[0].pT(), weight);
      _h["bjet_2_pT"]->fill(bjets[1].pT(), weight);
      _h["ljet_1_pT"]->fill(ljets[0].pT(), weight);
      _h["ljet_2_pT"]->fill(ljets[1].pT(), weight);

      // Construct the hadronically decaying W momentum 4-vector from pairs of
      // non-b-tagged jets. The pair which best matches the W mass is used. We start
      // with an always terrible 4-vector estimate which should always be "beaten" by
      // a real jet pair.
      FourMomentum W(10*(sqrtS()>0.?sqrtS():14000.), 0, 0, 0);
      for (size_t i = 0; i < ljets.size()-1; ++i) {
        for (size_t j = i + 1; j < ljets.size(); ++j) {
          const FourMomentum Wcand = ljets[i].momentum() + ljets[j].momentum();
          MSG_TRACE(i << "," << j << ": candidate W mass = " << Wcand.mass()/GeV
                    << " GeV, vs. incumbent candidate with " << W.mass()/GeV << " GeV");
          if (fabs(Wcand.mass() - 80.4*GeV) < fabs(W.mass() - 80.4*GeV)) {
            W = Wcand;
          }
        }
      }
      MSG_DEBUG("Candidate W mass = " << W.mass() << " GeV");

      // There are two b-jets with which this can be combined to make the
      // hadronically decaying top, one of which is correct and the other is
      // not... but we have no way to identify which is which, so we construct
      // both possible top momenta and fill the histograms with both.
      const FourMomentum t1 = W + bjets[0].momentum();
      const FourMomentum t2 = W + bjets[1].momentum();
      _h["W_mass"]->fill(W.mass(), weight);
      _h["t_mass"]->fill(t1.mass(), weight);
      _h["t_mass"]->fill(t2.mass(), weight);

      // Placing a cut on the well-known W mass helps to reduce backgrounds
      if (inRange(W.mass()/GeV, 75.0, 85.0)) {
        MSG_DEBUG("W found with mass " << W.mass()/GeV << " GeV");
        _h["t_mass_W_cut"]->fill(t1.mass(), weight);
        _h["t_mass_W_cut"]->fill(t2.mass(), weight);

        _h["jetb_1_jetb_2_dR"]->fill(deltaR(bjets[0].momentum(), bjets[1].momentum()),weight);
        _h["jetb_1_jetb_2_deta"]->fill(fabs(bjets[0].eta()-bjets[1].eta()),weight);
        _h["jetb_1_jetb_2_dphi"]->fill(deltaPhi(bjets[0].momentum(),bjets[1].momentum()),weight);

        _h["jetb_1_jetl_1_dR"]->fill(deltaR(bjets[0].momentum(), ljets[0].momentum()),weight);
        _h["jetb_1_jetl_1_deta"]->fill(fabs(bjets[0].eta()-ljets[0].eta()),weight);
        _h["jetb_1_jetl_1_dphi"]->fill(deltaPhi(bjets[0].momentum(),ljets[0].momentum()),weight);

        _h["jetl_1_jetl_2_dR"]->fill(deltaR(ljets[0].momentum(), ljets[1].momentum()),weight);
        _h["jetl_1_jetl_2_deta"]->fill(fabs(ljets[0].eta()-ljets[1].eta()),weight);
        _h["jetl_1_jetl_2_dphi"]->fill(deltaPhi(ljets[0].momentum(),ljets[1].momentum()),weight);

        _h["jetb_1_W_dR"]->fill(deltaR(bjets[0].momentum(), W),weight);
        _h["jetb_1_W_deta"]->fill(fabs(bjets[0].eta()-W.eta()),weight);
        _h["jetb_1_W_dphi"]->fill(deltaPhi(bjets[0].momentum(),W),weight);

        if (_mode > 0) {
          FourMomentum l=lfs.chargedLeptons()[0].momentum();
          _h["jetb_1_l_dR"]->fill(deltaR(bjets[0].momentum(), l),weight);
          _h["jetb_1_l_deta"]->fill(fabs(bjets[0].eta()-l.eta()),weight);
          _h["jetb_1_l_dphi"]->fill(deltaPhi(bjets[0].momentum(),l),weight);
          _h["jetb_1_l_mass"]->fill(FourMomentum(bjets[0].momentum()+l).mass(), weight);

          if (nLeps > 1) {
            FourMomentum l=lfs.chargedLeptons()[1].momentum();
            _h["jetb_1_l2_dR"]->fill(deltaR(bjets[0].momentum(), l),weight);
            _h["jetb_1_l2_deta"]->fill(fabs(bjets[0].eta()-l.eta()),weight);
            _h["jetb_1_l2_dphi"]->fill(deltaPhi(bjets[0].momentum(),l),weight);
            _h["jetb_1_l2_mass"]->fill(FourMomentum(bjets[0].momentum()+l).mass(), weight);
          }
        }
      }

    }


    void finalize() {
      const double sf = crossSection() / sumOfWeights();
      for (auto hist : _h) { scale(hist.second, sf); }
    }

    //@}

  protected:

      size_t _mode;


  private:

    // @name Histogram data members
    //@{
    map<string, Histo1DPtr> _h;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTBAR);
}
