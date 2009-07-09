// -*- C++ -*-
#ifndef RIVET_OPAL_1998_S3780481_HH
#define RIVET_OPAL_1998_S3780481_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief OPAL flavour dependent fragmentation paper
  /// @author Hendrik Hoeth
  class OPAL_1998_S3780481 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    OPAL_1998_S3780481();

    /// Factory method.
    static Analysis* create() { 
      return new OPAL_1998_S3780481(); 
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
    double _weightedTotalPartNum;

    double _SumOfudsWeights;
    double _SumOfcWeights;
    double _SumOfbWeights;

    AIDA::IHistogram1D *_histXpuds;
    AIDA::IHistogram1D *_histXpc;
    AIDA::IHistogram1D *_histXpb;
    AIDA::IHistogram1D *_histXpall;
    AIDA::IHistogram1D *_histLogXpuds;
    AIDA::IHistogram1D *_histLogXpc;
    AIDA::IHistogram1D *_histLogXpb;
    AIDA::IHistogram1D *_histLogXpall;
    AIDA::IHistogram1D *_histMultiChargeduds;
    AIDA::IHistogram1D *_histMultiChargedc;
    AIDA::IHistogram1D *_histMultiChargedb;
    AIDA::IHistogram1D *_histMultiChargedall;

    //@}

  };

}

#endif
