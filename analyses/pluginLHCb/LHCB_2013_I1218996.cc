// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// LHCb prompt charm hadron pT and rapidity spectra
  class LHCB_2013_I1218996 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    LHCB_2013_I1218996()
      : Analysis("LHCB_2013_I1218996")
    {    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections
      declare(UnstableParticles(), "UFS");

      /// Book histograms
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add( 2.0, 2.5, book(tmp, 3, 1, 1) );}
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add( 2.5, 3.0, book(tmp, 3, 1, 2) );}
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add( 3.0, 3.5, book(tmp, 3, 1, 3) );}
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add( 3.5, 4.0, book(tmp, 3, 1, 4) );}
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add( 4.0, 4.5, book(tmp, 3, 1, 5) );}

      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add( 2.0, 2.5, book(tmp, 2, 1, 1) );}
      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add( 2.5, 3.0, book(tmp, 2, 1, 2) );}
      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add( 3.0, 3.5, book(tmp, 2, 1, 3) );}
      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add( 3.5, 4.0, book(tmp, 2, 1, 4) );}
      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add( 4.0, 4.5, book(tmp, 2, 1, 5) );}

      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add( 2.0, 2.5, book(tmp, 5, 1, 1) );}
      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add( 2.5, 3.0, book(tmp, 5, 1, 2) );}
      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add( 3.0, 3.5, book(tmp, 5, 1, 3) );}
      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add( 3.5, 4.0, book(tmp, 5, 1, 4) );}
      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add( 4.0, 4.5, book(tmp, 5, 1, 5) );}

      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add( 2.0, 2.5, book(tmp, 4, 1, 1) );}
      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add( 2.5, 3.0, book(tmp, 4, 1, 2) );}
      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add( 3.0, 3.5, book(tmp, 4, 1, 3) );}
      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add( 3.5, 4.0, book(tmp, 4, 1, 4) );}
      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add( 4.0, 4.5, book(tmp, 4, 1, 5) );}

      book(_h_pdg4122_Lambdac_pT ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      /// @todo Use PrimaryHadrons to avoid double counting and automatically remove the contributions from unstable?
      const UnstableParticles &ufs = apply<UnstableParticles> (event, "UFS");
      for (const Particle& p : ufs.particles() ) {

        // We're only interested in charm hadrons
        if (!p.isHadron() || !p.hasCharm()) continue;

        // Kinematic acceptance
        const double y = p.absrap(); ///< Double analysis efficiency with a "two-sided LHCb"
        const double pT = p.pT();

        // Fiducial acceptance of the measurements
        if (pT > 8.0*GeV || y < 2.0 || y > 4.5) continue;

        /// Experimental selection removes non-prompt charm hadrons: we ignore those from b decays
        if (p.fromBottom()) continue;

        switch (p.abspid()) {
        case 411:
          _h_pdg411_Dplus_pT_y.fill(y, pT/GeV, weight);
          break;
        case 421:
          _h_pdg421_Dzero_pT_y.fill(y, pT/GeV, weight);
          break;
        case 431:
          _h_pdg431_Dsplus_pT_y.fill(y, pT/GeV, weight);
          break;
        case 413:
          _h_pdg413_Dstarplus_pT_y.fill(y, pT/GeV, weight);
          break;
        case 4122:
          _h_pdg4122_Lambdac_pT->fill(pT/GeV, weight);
          break;
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double scale_factor = 0.5 * crossSection()/microbarn / sumOfWeights();
      /// Avoid the implicit division by the bin width in the BinnedHistogram::scale method.
      for (Histo1DPtr h : _h_pdg411_Dplus_pT_y.histos()) h->scaleW(scale_factor);
      for (Histo1DPtr h : _h_pdg421_Dzero_pT_y.histos()) h->scaleW(scale_factor);
      for (Histo1DPtr h : _h_pdg431_Dsplus_pT_y.histos()) h->scaleW(scale_factor);
      for (Histo1DPtr h : _h_pdg413_Dstarplus_pT_y.histos()) h->scaleW(scale_factor);
      _h_pdg4122_Lambdac_pT->scaleW(scale_factor);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    BinnedHistogram _h_pdg411_Dplus_pT_y;
    BinnedHistogram _h_pdg421_Dzero_pT_y;
    BinnedHistogram _h_pdg431_Dsplus_pT_y;
    BinnedHistogram _h_pdg413_Dstarplus_pT_y;
    Histo1DPtr _h_pdg4122_Lambdac_pT;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(LHCB_2013_I1218996);

}
