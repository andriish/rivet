// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Measurement of D0 Run II Z \f$ p_\perp \f$ diff cross-section shape
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  /// @author Frank Siegert
  class D0_2007_S7075677 : public Analysis {

  public:

    /// Default constructor.
    D0_2007_S7075677() : Analysis("D0_2007_S7075677")
    {
      // Run II Z rapidity
    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      ZFinder zfinder(fs, -MAXRAPIDITY, MAXRAPIDITY, 0.0*GeV, PID::ELECTRON,
                      71.0*GeV, 111.0*GeV, 0.2, true, true);
      addProjection(zfinder, "ZFinder");

      _h_yZ = bookHisto1D(1, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = e.weight();

      const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() == 1) {
        const Particles& el(zfinder.constituents());
        if (el[0].momentum().pT() > 25.0*GeV || el[1].momentum().pT() > 25.0*GeV) {
          double yZ = fabs(zfinder.bosons()[0].momentum().rapidity());
          _h_yZ->fill(yZ, weight);
        }
      }
      else {
        MSG_DEBUG("No unique lepton pair found.");
      }
    }


    // Finalize
    void finalize() {
      // Data seems to have been normalized for the avg of the two sides
      // (+ve & -ve rapidity) rather than the sum, hence the 0.5:
      normalize(_h_yZ, 0.5);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_yZ;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2007_S7075677);

}
