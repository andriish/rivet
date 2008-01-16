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
    ExampleAnalysis()
      : _cmultproj(_cfsproj), _cnmultproj(_fsproj), 
        _thrustproj(_cfsproj), _sphericityproj(_cfsproj)
    {
      addProjection(_fsproj);
      addProjection(_cfsproj);
      addProjection(_cmultproj);
      addProjection(_cnmultproj);
      addProjection(_thrustproj);
      addProjection(_thrustproj);
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

    /// The final state projectors.
    FinalState _fsproj;
    ChargedFinalState _cfsproj;

    /// The multiplicity projectors (charged and all).
    Multiplicity _cmultproj;
    Multiplicity _cnmultproj;

    /// The thrust projector.
    Thrust _thrustproj;

    /// The sphericity projector.
    Sphericity _sphericityproj;


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
