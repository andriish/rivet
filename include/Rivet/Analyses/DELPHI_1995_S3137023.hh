// -*- C++ -*-
#ifndef RIVET_DELPHI_1995_S3137023_HH
#define RIVET_DELPHI_1995_S3137023_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// Implementation of DELPHI strange baryon paper
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
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "3137023";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Strange baryon production in Z hadronic decays at Delphi";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "DELPHI";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "1995";
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /// Hide the assignment operator
    DELPHI_1995_S3137023& operator=(const DELPHI_1995_S3137023&);

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
