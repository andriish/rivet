// -*- C++ -*-
#ifndef RIVET_DELPHI_1995_S3137023_HH
#define RIVET_DELPHI_1995_S3137023_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief DELPHI strange baryon paper
  /// @author Hendrik Hoeth
  class DELPHI_1995_S3137023 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    DELPHI_1995_S3137023();

    /// Factory method.
    static Analysis* create() { 
      return new DELPHI_1995_S3137023(); 
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
    double _weightedTotalNumXiMinus;
    double _weightedTotalNumSigma1385Plus;

    AIDA::IHistogram1D *_histXpXiMinus;
    AIDA::IHistogram1D *_histXpSigma1385Plus;
    //@}

  };

}

#endif
