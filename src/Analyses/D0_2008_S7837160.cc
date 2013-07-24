// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// @brief D0 Run II measurement of W charge asymmetry
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7837160 : public Analysis {
  public:

    /// Default constructor.
    D0_2008_S7837160()
      : Analysis("D0_2008_S7837160")
    {
      // Run II W charge asymmetry
    }


    /// @name Analysis methods
    //@{

    // Book histograms and set up projections
    void init() {
      // Projections
      /// @todo Use separate pT and ETmiss cuts in WFinder
      FinalState fs;
      const WFinder wfe(fs, -5, 5, 25.0*GeV, PID::ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      addProjection(wfe, "WFe");

      // Cross-section histograms
      _h_dsigplus_deta_25_35  = bookHisto1D(1,1,1,"/dsigplus_deta_25_35");
      _h_dsigminus_deta_25_35 = bookHisto1D(1,1,1,"/dsigminus_deta_25_35");
      _h_dsigplus_deta_35     = bookHisto1D(1,1,1,"/dsigplus_deta_35");
      _h_dsigminus_deta_35    = bookHisto1D(1,1,1,"/dsigminus_deta_35");
      _h_dsigplus_deta_25     = bookHisto1D(1,1,1,"/dsigplus_deta_25");
      _h_dsigminus_deta_25    = bookHisto1D(1,1,1,"/dsigminus_deta_25");

      _h_asym1 = bookScatter2D(1, 1, 1);
      _h_asym2 = bookScatter2D(1, 1, 2);
      _h_asym3 = bookScatter2D(1, 1, 3);
    }


    /// Do the analysis
    void analyze(const Event & event) {
      const WFinder& wf = applyProjection<WFinder>(event, "WFe");
      if (wf.bosons().size() == 0) {
        MSG_DEBUG("No W candidates found: vetoing");
        vetoEvent;
      }

      // Require that leptons have Et >= 25 GeV
      /// @todo Use pT cut in WFinder
      /// @todo Any ETmiss cut?
      FourMomentum p_e=wf.constituentLeptons()[0].momentum();
      int chg_e = PID::threeCharge(wf.constituentLeptons()[0].pdgId());
      if (p_e.eta() < 0) chg_e *= -1;
      assert(chg_e != 0);

      const double weight = event.weight();
      const double eta_e = fabs(p_e.eta());
      const double et_e = p_e.Et();
      if (et_e < 35*GeV) {
        // 25 <= ET < 35
        if (chg_e < 0) {
          _h_dsigminus_deta_25_35->fill(eta_e, weight);
        } else {
          _h_dsigplus_deta_25_35->fill(eta_e, weight);
        }
      } else {
        // ET >= 35
        if (chg_e < 0) {
          _h_dsigminus_deta_35->fill(eta_e, weight);
        } else {
          _h_dsigplus_deta_35->fill(eta_e, weight);
        }
      }
      // Inclusive: ET >= 25
      if (chg_e < 0) {
        _h_dsigminus_deta_25->fill(eta_e, weight);
      } else {
        _h_dsigplus_deta_25->fill(eta_e, weight);
      }
    }


    /// Finalize
    void finalize() {

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta
      // + dsig-/deta) for each Et region
      divide(*_h_dsigplus_deta_25_35 - *_h_dsigminus_deta_25_35,
             *_h_dsigplus_deta_25_35 + *_h_dsigminus_deta_25_35,
             _h_asym1);

      divide(*_h_dsigplus_deta_35 - *_h_dsigminus_deta_35,
             *_h_dsigplus_deta_35 + *_h_dsigminus_deta_35,
             _h_asym2);

      divide(*_h_dsigplus_deta_25 - *_h_dsigminus_deta_25,
             *_h_dsigplus_deta_25 + *_h_dsigminus_deta_25,
             _h_asym3);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_dsigplus_deta_25_35, _h_dsigminus_deta_25_35;
    Histo1DPtr _h_dsigplus_deta_35, _h_dsigminus_deta_35;
    Histo1DPtr _h_dsigplus_deta_25, _h_dsigminus_deta_25;


    Scatter2DPtr _h_asym1;
    Scatter2DPtr _h_asym2;
    Scatter2DPtr _h_asym3;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2008_S7837160);

}
