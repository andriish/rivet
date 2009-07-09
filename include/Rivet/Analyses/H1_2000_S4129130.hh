// -*- C++ -*-
#ifndef RIVET_H1_2000_S4129130_HH
#define RIVET_H1_2000_S4129130_HH

#include "Rivet/Analysis.hh" 

namespace Rivet {


  /// @brief H1 energy flow and charged particle spectra
  /// @author Peter Richardson
  /// Based on the HZtool analysis
  class H1_2000_S4129130 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    H1_2000_S4129130();

    /// Factory method.
    static Analysis* create() { 
      return new H1_2000_S4129130(); 
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
    vector<AIDA::IHistogram1D *> _histETLowQa;
    vector<AIDA::IHistogram1D *> _histETHighQa;
    vector<AIDA::IHistogram1D *> _histETLowQb;
    vector<AIDA::IHistogram1D *> _histETHighQb;
    AIDA::IProfile1D * _histAverETCentral;
    AIDA::IProfile1D * _histAverETFrag;
    //@}

    /// @name storage of weights for normalisation
    //@{
    vector<double> _weightETLowQa;
    vector<double> _weightETHighQa;
    vector<double> _weightETLowQb;
    vector<double> _weightETHighQb;
    //@}
  };

}

#endif
