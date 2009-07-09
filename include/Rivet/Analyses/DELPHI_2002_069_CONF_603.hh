// -*- C++ -*-
#ifndef RIVET_DELPHI_2002_069_CONF_603_HH
#define RIVET_DELPHI_2002_069_CONF_603_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief DELPHI b-fragmentation measurement
  /// @author Hendrik Hoeth
  class DELPHI_2002_069_CONF_603 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    DELPHI_2002_069_CONF_603();

    /// Factory method.
    static Analysis* create() { 
      return new DELPHI_2002_069_CONF_603(); 
    }

    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the 
    /// inclusive single particle distributions' normalisations.

    AIDA::IHistogram1D *_histXbprim;
    AIDA::IHistogram1D *_histXbweak;

    AIDA::IProfile1D *_histMeanXbprim;
    AIDA::IProfile1D *_histMeanXbweak;

    //@}

  };

}

#endif
