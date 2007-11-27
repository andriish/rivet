// -*- C++ -*-
#ifndef RIVET_ExampleAnalysis_HH
#define RIVET_ExampleAnalysis_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This class just measures the charged multiplicity
  class ExampleAnalysis : public Analysis {

  public:

    /// Default constructor
    ExampleAnalysis()
      : _multproj(_fsproj), _thrustproj(_fsproj)
    {
      addProjection(_fsproj);
      addProjection(_multproj);
      addProjection(_thrustproj);
    }

    /// Factory method
    static Analysis* create() { 
      return new ExampleAnalysis(); 
    }

    /// Get the name of this analysis.
    string getName() const {
      return "Test";
    }

    void init();
    
    void analyze(const Event& event);
    
    void finalize();

  private:

    /// The final state projector.
    FinalState _fsproj;

    /// The multiplicity projector.
    Multiplicity _multproj;

    /// The thrust projector.
    Thrust _thrustproj;

  private:

    /// Hide the assignment operator
    ExampleAnalysis& operator=(const ExampleAnalysis&);

    //@{
    /// Histograms
    AIDA::IHistogram1D* _histTot;
    AIDA::IHistogram1D* _histChTot;
    AIDA::IHistogram1D* _histUnchTot;
    AIDA::IHistogram1D* _histHadrTot;
    AIDA::IHistogram1D* _histHadrChTot;
    AIDA::IHistogram1D* _histHadrUnchTot;
    AIDA::IHistogram1D* _histThrust;
    //@}

  };

}

#endif
