// -*- C++ -*-
#ifndef RIVET_D0_2009_S8320160_HH
#define RIVET_D0_2009_S8320160_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class D0_2009_S8320160 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    D0_2009_S8320160();

    /// Factory method 
    static Analysis* create() {
      return new D0_2009_S8320160();
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
    BinnedHistogram<double> _h_chi_dijet;
    //@}
    
  };


}

#endif
