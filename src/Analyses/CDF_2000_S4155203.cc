// -*- C++ -*-
// CDF Z pT analysis

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /*
   * @brief CDF Run I Z pT in Drell-Yan events
   * @author Hendrik Hoeth
   */ 
  class CDF_2000_S4155203 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV.
    CDF_2000_S4155203() 
      : Analysis("CDF_2000_S4155203")
    {
      setBeams(PROTON, ANTIPROTON);
      ZFinder zfinder(FinalState(), ELECTRON, 66.0*GeV, 116.0*GeV, 0.2);
      addProjection(zfinder, "ZFinder");
    }
    
    //@}


    /// @name Analysis methods
    //@{
    
    /// Book histograms
    void init() {
      _hist_zpt = bookHistogram1D(1, 1, 1);
    }
    
    
    /// Do the analysis
    void analyze(const Event& e) {
      const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
      if (zfinder.particles().size() != 1) {
        getLog() << Log::DEBUG << "No unique e+e- pair found" << endl;
        vetoEvent;
      }
      
      FourMomentum pZ = zfinder.particles()[0].momentum();    
      getLog() << Log::DEBUG << "Dilepton mass = " << pZ.mass()/GeV << " GeV"  << endl;
      getLog() << Log::DEBUG << "Dilepton pT   = " << pZ.pT()/GeV << " GeV" << endl;
      _hist_zpt->fill(pZ.pT()/GeV, e.weight());
    }
    
    
    void finalize() {
      // Normalize to the experimental cross-section
      /// @todo Get norm from generator cross-section
      normalize(_hist_zpt, 247.4);
    }
    
    //@}


  private:

    // Histogram
    AIDA::IHistogram1D *_hist_zpt;

  };

    

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2000_S4155203> plugin_CDF_2000_S4155203;

}
