// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"

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
    	// select phi mesons in fiducial phase-space (symmetrical LHCb)
      declare(UnstableParticles(Cuts::abspid == 333 && Cuts::absrapIn(2.44, 4.06) && Cuts::ptIn(0.6*GeV, 5.0*GeV)), "phiFS");

      {Histo1DPtr tmp; _h_Phi_pT_y.add(  2.44, 2.62, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_Phi_pT_y.add(  2.62, 2.80, book(tmp, 2, 1, 2));}
      {Histo1DPtr tmp; _h_Phi_pT_y.add(  2.80, 2.98, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_Phi_pT_y.add(  2.98, 3.16, book(tmp, 3, 1, 2));}
      {Histo1DPtr tmp; _h_Phi_pT_y.add(  3.16, 3.34, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_Phi_pT_y.add(  3.34, 3.52, book(tmp, 4, 1, 2));}
      {Histo1DPtr tmp; _h_Phi_pT_y.add(  3.52, 3.70, book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_Phi_pT_y.add(  3.70, 3.88, book(tmp, 5, 1, 2));}
      {Histo1DPtr tmp; _h_Phi_pT_y.add(  3.88, 4.06, book(tmp, 6, 1, 1));}
      book(_h_Phi_pT ,7, 1, 1);
      book(_h_Phi_y ,8, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze (const Event& event) {
      const double weight = 1.;
      const UnstableParticles& ufs = apply<UnstableParticles> (event, "phiFS");

      for (const Particle& p : ufs.particles()) {
      	double y  = p.rapidity();
        double pT = p.pt();
        _h_Phi_y->fill (y, weight);
        _h_Phi_pT->fill (pT/MeV, weight);
        _h_Phi_pT_y.fill(y, pT/GeV, weight);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
    	// correct for symmetrical detector
      double scale_factor = crossSectionPerEvent()/microbarn/2.;
      scale (_h_Phi_y, scale_factor);
      scale (_h_Phi_pT, scale_factor);
      _h_Phi_pT_y.scale(scale_factor/1000., this);
    }

    //@}

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Phi_y;
    Histo1DPtr _h_Phi_pT;
    BinnedHistogram _h_Phi_pT_y;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(LHCB_2011_I919315);

}

//@}
