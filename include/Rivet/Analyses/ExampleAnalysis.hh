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

  /// This class just measures a few random things as an example.
  class ExampleAnalysis : public Analysis {

  public:

    /// Default constructor
    ExampleAnalysis() {
      const FinalState& cnfs = addProjection(*new FinalState(), "FS");
      const ChargedFinalState& cfs = addProjection(*new ChargedFinalState(), "CFS");
      addProjection(*new Multiplicity(cfs), "CMult");
      addProjection(*new Multiplicity(cnfs), "CNMult");
      addProjection(*new Thrust(cfs), "Thrust");
      addProjection(*new Sphericity(cfs), "Sphericity");
    }

    /// Factory method
    static Analysis* create() { 
      return new ExampleAnalysis(); 
    }

  public:

    /// @name Publication metadata
    //@{
    /// Return the name of this analysis
    string getName() const {
      return "Example";
    }
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "NONE";
    }
    /// Get a description of the analysis.
    // string getDescription() const {
    //   return "";
    // }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "NONE";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
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

    /// Hide the assignment operator
    ExampleAnalysis& operator=(const ExampleAnalysis&);

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
