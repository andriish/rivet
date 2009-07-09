// -*- C++ -*-
#ifndef RIVET_D0_2009_S8202443_HH
#define RIVET_D0_2009_S8202443_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class D0_2009_S8202443 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    D0_2009_S8202443();

    /// Factory method 
    static Analysis* create() {
      return new D0_2009_S8202443();
    }
    //@}


    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_jet1_pT;
    AIDA::IHistogram1D * _h_jet2_pT;
    AIDA::IHistogram1D * _h_jet3_pT;
    AIDA::IHistogram1D * _h_jet1_pT_constrained;
    AIDA::IHistogram1D * _h_jet2_pT_constrained;
    AIDA::IHistogram1D * _h_jet3_pT_constrained;
    //@}
    
    double _sum_of_weights, _sum_of_weights_constrained;

  };


}

#endif
