// -*- C++ -*-
#ifndef RIVET_CDF_2008_S8093652_HH
#define RIVET_CDF_2008_S8093652_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CDF_2008_S8093652 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    CDF_2008_S8093652();

    /// Factory method 
    static Analysis* create() {
      return new CDF_2008_S8093652();
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
    AIDA::IHistogram1D* _h_m_dijet;
    //@}
    
  };


}

#endif
