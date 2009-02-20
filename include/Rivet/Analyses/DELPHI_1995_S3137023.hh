// -*- C++ -*-
#ifndef RIVET_DELPHI_1995_S3137023_HH
#define RIVET_DELPHI_1995_S3137023_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief DELPHI strange baryon paper
  /// @author Hendrik Hoeth
  class DELPHI_1995_S3137023 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    DELPHI_1995_S3137023() 
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");
      addProjection(UnstableFinalState(), "UFS");
      _weightedTotalNumXiMinus = 0;
      _weightedTotalNumSigma1385Plus = 0;
    }


    /// Factory method.
    static Analysis* create() { 
      return new DELPHI_1995_S3137023(); 
    }

    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "3137023";
    }
    /// A short description of the analysis.
    string summary() const {
      return "Strange baryon production in Z hadronic decays at Delphi";
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Measurement of the $\\Xi^-$ and $\\Sigma^+(1385)/\\Sigma^-(1385)$ scaled momentum "
         << "distributions by DELPHI at LEP 1. The paper also has the production "
	 << "cross-sections of these particles, but that's not implemented in Rivet.";
      return os.str();
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "DELPHI";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "LEP 1";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1995";
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
      ret += "Z. Phys. C, 67, 543--554 (1995)";
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
    double _weightedTotalNumXiMinus;
    double _weightedTotalNumSigma1385Plus;

    AIDA::IHistogram1D *_histXpXiMinus;
    AIDA::IHistogram1D *_histXpSigma1385Plus;
    //@}

  };

}

#endif
