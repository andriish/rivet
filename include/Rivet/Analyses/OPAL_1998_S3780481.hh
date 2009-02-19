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
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurement of scaled momentum distributions and total "
         << "charged multiplicities in flavour tagged events at LEP 1. "
         << "OPAL measured these observables in uds-, c-, and b-events "
         << "separately. An inclusive measurement is also included.";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "OPAL";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "LEP1";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1998";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Hendrik Hoeth <hendrik.hoeth@cern.ch>";
      return ret;
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "Hadronic Z decay events generated on the Z pole (sqrt(s) = 91.2 GeV)";
      return os.str();
    }
    string status() const {
      return "VALIDATED";
    }
    /// Journal, and preprint references
    vector<string> references() const {
      vector<string> ret;
      ret += "Eur. Phys. J, C7, 369--381 (1999)"; 
      ret += "hep-ex/9807004";
      return ret;
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
