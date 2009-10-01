// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  class CDF_2009_S8383952 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2009_S8383952()
      : Analysis("CDF_2009_S8383952") 
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections here
      ZFinder zfinder(-2.8, 2.8, 20.0*GeV, ELECTRON,
                      66.0*GeV, 116.0*GeV, 0.2);
      addProjection(zfinder, "ZFinder");


      /// Book histograms here, e.g.:
      _h_yZ = bookHistogram1D(1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
      if (zfinder.particles().size() == 1) {
        double yZ = fabs(zfinder.particles()[0].momentum().rapidity());
        _h_yZ->fill(yZ, weight);
      }
      else {
        getLog() << Log::DEBUG << "no unique lepton pair found." << endl;
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_yZ, crossSection()/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{

    AIDA::IHistogram1D *_h_yZ;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2009_S8383952> plugin_CDF_2009_S8383952;


}
