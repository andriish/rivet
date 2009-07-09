// -*- C++ -*-
#ifndef RIVET_ALEPH_1991_S2435284_HH
#define RIVET_ALEPH_1991_S2435284_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of ALEPH LEP1 charged multiplicity
  /// @author Andy Buckley
  class ALEPH_1991_S2435284 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor.
    ALEPH_1991_S2435284();

    /// Factory method.
    static Analysis* create() { 
      return new ALEPH_1991_S2435284(); 
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event & event);
    virtual void finalize();
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histChTot;
    //@}

  };

}


#endif
