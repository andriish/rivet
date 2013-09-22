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
    {    }


    /// @name Analysis methods
    //@{

    // Book histograms and set up projections
    void init() {
      // Projections
      FinalState fs;
      /// @todo Use separate pT and ETmiss cuts in WFinder
      const WFinder wfe(fs, -5, 5, 25.0*GeV, PID::ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      addProjection(wfe, "WFe");

      // Cross-section histograms
      for (size_t pmindex = 0; pmindex <= 1; ++pmindex) {
        _hs_dsigpm_deta_25_35[pmindex] = bookHisto1D(1, 1, 1, "/TMP/dsigpm_deta_25_35_" + to_str(pmindex));
        _hs_dsigpm_deta_35[pmindex] = bookHisto1D(1, 1, 2, "/TMP/dsigpm_deta_35_" + to_str(pmindex));
        _hs_dsigpm_deta_25[pmindex] = bookHisto1D(1, 1, 3, "/TMP/dsigpm_deta_25_" + to_str(pmindex));
      }

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
      FourMomentum p_e = wf.constituentLeptons()[0].momentum();
      int chg_e = PID::threeCharge(wf.constituentLeptons()[0].pdgId());
      if (p_e.eta() < 0) chg_e *= -1;
      assert(chg_e != 0);

      const double weight = event.weight();
      const double eta_e = fabs(p_e.eta());
      const double et_e = p_e.Et();

      // Fill histos with appropriate +- indexing
      const size_t pmindex = (chg_e > 0) ? 0 : 1;
      if (et_e < 35*GeV) _hs_dsigpm_deta_25_35[pmindex]->fill(eta_e, weight);
      else _hs_dsigpm_deta_35[pmindex]->fill(eta_e, weight);
      _hs_dsigpm_deta_25[pmindex]->fill(eta_e, weight);
    }


    /// @name Helper functions for constructing asymmetry histograms in finalize()
    //@{
    void calc_asymm(const Histo1DPtr plus, const Histo1DPtr minus, Scatter2DPtr target) {
      divide(*plus - *minus, *plus + *minus, target);
    }
    void calc_asymm(const Histo1DPtr histos[2], Scatter2DPtr target) {
      calc_asymm(histos[0], histos[1], target);
    }
    //@}


    /// @brief Finalize
    ///
    /// Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta) for each ET region
    void finalize() {
      calc_asymm(_hs_dsigpm_deta_25_35, _h_asym1);
      calc_asymm(_hs_dsigpm_deta_35, _h_asym2);
      calc_asymm(_hs_dsigpm_deta_25, _h_asym3);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hs_dsigpm_deta_25_35[2], _hs_dsigpm_deta_35[2], _hs_dsigpm_deta_25[2];
    Scatter2DPtr _h_asym1, _h_asym2, _h_asym3;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2008_S7837160);

}
