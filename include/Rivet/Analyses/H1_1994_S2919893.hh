// -*- C++ -*-
#ifndef RIVET_H1_1994_S2919893_HH
#define RIVET_H1_1994_S2919893_HH

#include "Rivet/Analysis.hh" 

namespace Rivet {


  /// @brief H1 energy flow and charged particle spectra
  /// @author Peter Richardson
  /// Based on the equivalent HZTool analysis
  class H1_1994_S2919893 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    H1_1994_S2919893();

    /// Factory method.
    static Analysis* create() { 
      return new H1_1994_S2919893(); 
    }

    //@}
    

    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /**
     *  Polar angle with right direction of the beam
     */
    inline double beamAngle(const FourVector& v, const bool & order) {
      double thel = v.polarAngle()/degree;
      if(thel<0.) thel+=180.;
      if(!order) thel = 180.-thel;
      return thel;
    }

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histEnergyFlowLowX;
    AIDA::IHistogram1D *_histEnergyFlowHighX;
    AIDA::IHistogram1D *_histEECLowX;
    AIDA::IHistogram1D *_histEECHighX;
    AIDA::IHistogram1D *_histSpectraW77;
    AIDA::IHistogram1D *_histSpectraW122;
    AIDA::IHistogram1D *_histSpectraW169;
    AIDA::IHistogram1D *_histSpectraW117;
    AIDA::IProfile1D *_histPT2;
    //@}

    /// @name storage of weight to calculate averages for normalisation
    //@{
    pair<double,double> _w77,_w122,_w169,_w117,_wEnergy;
    //@}
  };

}

#endif
