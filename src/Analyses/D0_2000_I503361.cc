// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  using namespace Cuts;

  /// @ D0 Run I Z \f$ p_\perp \f$ in Drell-Yan events
  /// @author Simone Amoroso
  class D0_2000_I503361 : public Analysis {
  public:

    /// Constructor
    D0_2000_I503361()
      : Analysis("D0_2000_I503361")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      ///  Initialise and register projections here
      ZFinder zfinder(FinalState(), Cuts::open(), PID::ELECTRON, 75*GeV, 105*GeV, 25.0*GeV, ZFinder::NOCLUSTER);
      addProjection(zfinder, "ZFinder");


      ///  Book histograms here, e.g.:
      //      _hist_z_xs = bookHisto1D(1, 1, 1);
      _hist_zpt = bookHisto1D(1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here
      const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1) {
        MSG_DEBUG("Num e+ e- pairs found = " << zfinder.bosons().size());
	vetoEvent;
      }
      const FourMomentum& pZ = zfinder.bosons()[0].momentum();
      if (pZ.mass2() < 0) {
	MSG_DEBUG("Negative Z mass**2 = " << pZ.mass2()/GeV2 << "!");
	vetoEvent;
      }

      MSG_DEBUG("Dilepton mass = " << pZ.mass()/GeV << " GeV");
      _hist_zpt->fill(pZ.pT(), weight);

      /// A posteriori trigger on both electrons to have E_T > 25 GeV
      if (zfinder.constituents()[0].Et()/GeV < 25 || 
          zfinder.constituents()[1].Et()/GeV < 25 ) {
        MSG_DEBUG("E_T trigger failed: l1 = " 
            <<  zfinder.constituents()[0].Et()/GeV << " GeV"
            << " l2 = " << zfinder.constituents()[1].Et()/GeV);
        vetoEvent;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_zpt, crossSection()/picobarn/sumOfWeights());
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


    /// @name Histograms
    //@{
    Histo1DPtr _hist_zpt;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2000_I503361);

}
