// -*- C++ -*-
#ifndef RIVET_ZEUS_2001_S4815815_HH
#define RIVET_ZEUS_2001_S4815815_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief ZEUS dijet photoproduction study used in the ZEUS Jets PDF fit
  ///
  /// This class is a reproduction of the HZTool routine for the ZEUS 
  /// dijet photoproduction paper which was used in the ZEUS Jets PDF fit.  
  ///
  /// @author Jon Butterworth
  class ZEUS_2001_S4815815 : public Analysis {

  public:

    /// Default constructor.
    ZEUS_2001_S4815815();

    /// Factory method.
    static Analysis* create() { 
      return new ZEUS_2001_S4815815(); 
    }

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetEt1;
    //@}


  };


}

#endif
