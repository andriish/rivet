// -*- C++ -*-
#ifndef RIVET_CDF_1988_S1865951_HH
#define RIVET_CDF_1988_S1865951_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    
  class CDF_1988_S1865951 : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// constructor
    CDF_1988_S1865951();

    /// Factory method
    static Analysis* create() { 
      return new CDF_1988_S1865951(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

    private:
    AIDA::IHistogram1D* _hist_pt630;
    AIDA::IHistogram1D* _hist_pt1800;
  };

}

#endif

