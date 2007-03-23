// -*- C++ -*-
#ifndef RIVET_TestAnalysis_H
#define RIVET_TestAnalysis_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This class just measures the charged multiplicity
  class TestAnalysis : public Analysis {

  public:

    /// Default constructor
    inline TestAnalysis() 
      : p_mult(p_fs), p_thrust(p_fs)
    { }

    /// Copy constructor.
    inline TestAnalysis(const TestAnalysis& x)
      : Analysis(x), p_fs(x.p_fs), p_mult(x.p_mult), p_thrust(x.p_thrust) 
    { }

    /// Destructor
    inline ~TestAnalysis() 
    { }

  public:

    /// The name of this analysis is "Test"
    inline std::string name() const {
      return "Test";
    }

  public:

    void init();
    
    void analyze(const Event & event);
    
    void finalize();

    /// Return the RivetInfo object of this analysis object.
    RivetInfo getInfo() const;

  private:

    /// The FinalState projector used by this analysis.
    FinalState p_fs;

    /// The Multiplicity projector used by this analysis.
    Multiplicity p_mult;

    /// The thrust projector
    Thrust p_thrust;

  private:

    /// Hide the assignment operator
    TestAnalysis & operator=(const TestAnalysis& x);

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
