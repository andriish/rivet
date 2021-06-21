// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Measurements of production cross sections of WZ and same-sign WW boson pairs in association with two jets in proton-proton collisions at 13 TeV
  class CMS_2020_I1794169 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2020_I1794169);

    float totalEvents = 0;
    float sumSelWWEvents = 0;
    float sumSelWZEvents = 0;
 
    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      _mode = 0;
      if ( getOption("LMODE") == "WZ" ) _mode = 1;

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);
      const FinalState fsjet4p7(Cuts::abseta < 4.7);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jet4p7fs(fsjet4p7, FastJets::ANTIKT, 0.4);
      declare(jet4p7fs, "jets4p7");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      bare_leps.acceptTauDecays(false);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      book(_h_WW_mjj   ,  9, 1, 1);
      book(_h_WW_mll   , 11, 1, 1);
      book(_h_WW_ptlmax, 13, 1, 1);
      book(_h_WZ_mjj   , 15, 1, 1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      totalEvents++;

      // Retrieve dressed leptons, sorted by pT
      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // Apply a lepton size requirement
      if (leptons.size() <= 1 || leptons.size() >= 4) return;

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets50 = apply<FastJets>(event, "jets4p7").jetsByPt(Cuts::pT > 50*GeV);

      // Remove all jets within dR < 0.4 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets50, leptons, 0.4);

      // Apply a njets >= 2 cut
      if (jets50.size() < 2) return;

      FourMomentum dijetCand = jets50[0].momentum() + jets50[1].momentum();
      double deltaEtaJJ = std::abs(jets50[0].eta()-jets50[1].eta());

      // Apply a mjj > 500 and detajj > 2.5 cuts
      if (dijetCand.mass() <= 500 || deltaEtaJJ <= 2.5) return;

      // W+W+ selection
      if (leptons.size() == 2 && leptons[0].pid() * leptons[1].pid() > 0 && _mode == 0) {
        FourMomentum dilCand = leptons[0].momentum() + leptons[1].momentum();
        if(dilCand.mass() > 20){
          double ptlmax = leptons[0].pt(); double ptlmin = leptons[1].pt();
          if(ptlmax < ptlmin) {
            ptlmax = leptons[1].pt(); ptlmin = leptons[0].pt();
          }

          _h_WW_mjj   ->fill(min(dijetCand.mass(),2999.999));
          _h_WW_mll   ->fill(min(dilCand.mass(),499.999));
          _h_WW_ptlmax->fill(min(ptlmax,299.999));

          sumSelWWEvents++;
        }
      }
      // WZ selection
      else if (leptons.size() == 3 && _mode == 1) {
        double mllZ = 10000; int iW = -1;
        if(leptons[0].pid() * leptons[1].pid() < 0 && std::abs(leptons[0].pid()) == std::abs(leptons[1].pid()) &&
           fabs((leptons[0].momentum() + leptons[1].momentum()).mass()-91.1876) < fabs(mllZ-91.1876)) {
          mllZ = (leptons[0].momentum() + leptons[1].momentum()).mass(); iW = 2;
        }

        if(leptons[0].pid() * leptons[2].pid() < 0 && std::abs(leptons[0].pid()) == std::abs(leptons[2].pid()) &&
           fabs((leptons[0].momentum() + leptons[2].momentum()).mass()-91.1876) < fabs(mllZ-91.1876)) {
          mllZ = (leptons[0].momentum() + leptons[2].momentum()).mass(); iW = 1;
        }

        if(leptons[1].pid() * leptons[2].pid() < 0 && std::abs(leptons[1].pid()) == std::abs(leptons[2].pid()) &&
           fabs((leptons[1].momentum() + leptons[2].momentum()).mass()-91.1876) < fabs(mllZ-91.1876)) {
          mllZ = (leptons[1].momentum() + leptons[2].momentum()).mass(); iW = 0;
        }

        if(iW >= 0 && fabs(mllZ-91.1876) < 15){
          _h_WZ_mjj->fill(min(dijetCand.mass(),2999.999));
 
          sumSelWZEvents++;
        }
      }

    }

    void normalizeToSum(Histo1DPtr hist) {
      double sum = 0.;
      for (size_t i = 0; i < hist->numBins(); ++i) {
        sum += hist->bin(i).height();
        float width = hist->bin(i).width();
        hist->bin(i).scaleW(width != 0 ? width : 1.);
      }
      if(hist->integral() > 0) scale(hist, 1./hist->integral());
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = (sumOfWeights() != 0) ? crossSection()/femtobarn/sumOfWeights() : 1.0;

      scale(_h_WW_mjj   , norm);
      scale(_h_WW_mll   , norm);
      scale(_h_WW_ptlmax, norm);
      scale(_h_WZ_mjj   , norm);

    }

    //@}

  protected:

    size_t _mode;

  private:
    /// @name Histograms
    //@{
    Histo1DPtr _h_WW_mjj, _h_WW_mll, _h_WW_ptlmax, _h_WZ_mjj;
    //@}


  };


  DECLARE_RIVET_PLUGIN(CMS_2020_I1794169);

}
