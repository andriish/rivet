// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  class LHCB_2011_I919315 : public Analysis {
  public:
    /// @name Constructors etc.
    //@{

    /// Constructor
    LHCB_2011_I919315()
      : Analysis("LHCB_2011_I919315")
    {
    }

    //@}
  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      addProjection(UnstableFinalState(), "UFS");

      _h_Phi_pT_y.addHistogram(  2.44, 2.62, bookHistogram1D(1, 1, 1));
      _h_Phi_pT_y.addHistogram(  2.62, 2.80, bookHistogram1D(1, 1, 2));
      _h_Phi_pT_y.addHistogram(  2.80, 2.98, bookHistogram1D(1, 1, 3));
      _h_Phi_pT_y.addHistogram(  2.98, 3.16, bookHistogram1D(1, 1, 4));
      _h_Phi_pT_y.addHistogram(  3.16, 3.34, bookHistogram1D(1, 1, 5));
      _h_Phi_pT_y.addHistogram(  3.34, 3.52, bookHistogram1D(1, 1, 6));
      _h_Phi_pT_y.addHistogram(  3.52, 3.70, bookHistogram1D(1, 1, 7));
      _h_Phi_pT_y.addHistogram(  3.70, 3.88, bookHistogram1D(1, 1, 8));
      _h_Phi_pT_y.addHistogram(  3.88, 4.06, bookHistogram1D(1, 1, 9));
      _h_Phi_pT = bookHistogram1D(1, 1, 10);
      _h_Phi_y = bookHistogram1D(1, 2, 1);
    }


    /// Perform the per-event analysis
    void analyze (const Event& event) {
      const double weight = event.weight();
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState> (event, "UFS");

      foreach (const Particle& p, ufs.particles()) {
        const PdgId id = abs(p.pdgId());

        if (id == 333) { // id 333 = phi-meson
          double y  = p.momentum().rapidity();
          double pT = p.momentum().perp();

          if (pT < 0.6*GeV || pT > 5.0*GeV || y < 2.44 || y > 4.06) {
            continue;
          }
          
          _h_Phi_y->fill (y, weight);
          _h_Phi_pT->fill (pT/MeV, weight);
          _h_Phi_pT_y.fill(y, pT/MeV, weight);
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      double scale_factor = crossSectionPerEvent()/microbarn;
      scale (_h_Phi_y, scale_factor);
      scale (_h_Phi_pT, scale_factor);
      _h_Phi_pT_y.scale(scale_factor, this);
    }

    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_Phi_y;
    AIDA::IHistogram1D * _h_Phi_pT;
    BinnedHistogram<double> _h_Phi_pT_y;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2011_I919315);

}

//@}
