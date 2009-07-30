// -*- C++ -*-
#ifndef RIVET_UA5_1982_S875503_HH
#define RIVET_UA5_1982_S875503_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    
  class UA5_1982_S875503 : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// constructor
    UA5_1982_S875503();

    /// Factory method
    static Analysis* create() { 
      return new UA5_1982_S875503(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

    private:
    
    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_etapp;
    AIDA::IHistogram1D* _hist_etappbar;
    //@}

  };

}

#endif

