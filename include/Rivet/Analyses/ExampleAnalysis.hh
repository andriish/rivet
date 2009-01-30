// -*- C++ -*-
#ifndef RIVET_ExampleAnalysis_HH
#define RIVET_ExampleAnalysis_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// @brief Just measures a few random things as an example.
  class ExampleAnalysis : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// Default constructor
    ExampleAnalysis();

    /// Factory method
    static Analysis* create() { 
      return new ExampleAnalysis(); 
    }

    //@}


  public:

    /// @name Publication metadata
    //@{
    /// Return the name of this analysis
    string name() const {
      return "Example";
    }
    /// Get a description of the analysis.
    string spiresId() const {
      return "NONE";
    }
    /// Get a description of the analysis.
    // string description() const {
    //   return "";
    // }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "NONE";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "NONE";
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:
    
    //@{
    /// Histograms
    AIDA::IHistogram1D* _histTot;
    AIDA::IHistogram1D* _histChTot;
    AIDA::IHistogram1D* _histHadrTot;
    AIDA::IHistogram1D* _histHadrChTot;
    AIDA::IHistogram1D* _histThrust;
    AIDA::IHistogram1D* _histMajor;
    AIDA::IHistogram1D* _histSphericity;
    AIDA::IHistogram1D* _histAplanarity;
    //@}

  };

}

#endif
