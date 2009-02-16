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
    /// A short description of the analysis.
    string spiresId() const {
      return "NONE";
    }
    /// A short description of the analysis.
    string summary() const {
      return "A demo to show aspects of writing a Rivet analysis";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "NONE";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      return "NONE";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "NONE";
    }
    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Andy Buckley <andy.buckley@durham.ac.uk>";
      return rtn;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "This analysis is a demonstration of the Rivet analysis structure "
         << "and functionality: booking histograms; the initialisation, analysis "
         << "and finalisation phases; and a simple loop over event particles. "
         << "It has no physical meaning, but can be used as a simple pedagogical "
         << "template for writing real analyses.";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "All event types will be accepted.";
      return os.str();
    }
    string status() const {
      return "EXAMPLE";
    }
    /// No journal or preprint references: this is a demo.
    vector<string> references() const {
      vector<string> ret;
      return ret;
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
