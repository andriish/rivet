// -*- C++ -*-
#ifndef RIVET_OPAL_1998_S3780481_HH
#define RIVET_OPAL_1998_S3780481_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief OPAL flavour dependent fragmentation paper
  /// @author Hendrik Hoeth
  class OPAL_1998_S3780481 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    OPAL_1998_S3780481() 
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");
      addProjection(InitialQuarks(), "IQF");
      _weightedTotalPartNum = 0;
      _SumOfudsWeights = 0;
      _SumOfcWeights = 0;
      _SumOfbWeights = 0;
    }


    /// Factory method.
    static Analysis* create() { 
      return new OPAL_1998_S3780481(); 
    }

    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "3780481";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Measurements of flavor dependent fragmentation functions in Z0 --> q anti-q events.";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "OPAL";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1998";
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
