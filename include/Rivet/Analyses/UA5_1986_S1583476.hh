// -*- C++ -*-
#ifndef RIVET_UA5_1986_S1583476_HH
#define RIVET_UA5_1986_S1583476_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    
  class UA5_1986_S1583476 : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// constructor
    UA5_1986_S1583476();

    /// Factory method
    static Analysis* create() { 
      return new UA5_1986_S1583476(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


    private:
    
    AIDA::IHistogram1D* _hist_eta200;
    AIDA::IHistogram1D* _hist_eta900;

  };

}

#endif

