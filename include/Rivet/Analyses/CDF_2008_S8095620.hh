// -*- C++ -*-
#ifndef RIVET_CDF_2008_S8095620_HH
#define RIVET_CDF_2008_S8095620_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// Implementation of CDF Run II Z+b-jet cross section paper
  class CDF_2008_S8095620 : public Analysis {

  public:


    /// @name Constructors etc.
    //@{

    /// Constructor.
    /// jet cuts: |eta| <= 1.5
    CDF_2008_S8095620();

    /// Factory method
    static Analysis* create() { 
      return new CDF_2008_S8095620(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:

    double _Rjet;
    double _JetPtCut;
    double _JetEtaCut;
    double _sumWeightSelected;
 
    //@{
    /// Histograms
    AIDA::IHistogram1D* _dSdET;
    AIDA::IHistogram1D* _dSdETA;
    AIDA::IHistogram1D* _dSdNJet; 
    AIDA::IHistogram1D* _dSdNbJet; 
    AIDA::IHistogram1D* _dSdZpT; 

    //@}

  };

}

#endif
