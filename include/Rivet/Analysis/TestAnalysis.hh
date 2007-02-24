// -*- C++ -*-
#ifndef RIVET_TestAnalysis_H
#define RIVET_TestAnalysis_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {

  /// This class just measures the charged multiplicity
  class TestAnalysis : public Analysis {

  public:

    /// Default constructor.
    inline TestAnalysis() : mult(fsproj) {}

    /// Copy constructor.
    inline TestAnalysis(const TestAnalysis& x)
      : Analysis(x), fsproj(x.fsproj), mult(x.mult) {}

    /// Destructor
    ~TestAnalysis();

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
    FinalState fsproj;

    /// The Multiplicity projector used by this analysis.
    Multiplicity mult;

  private:

    /// Hide the assignment operator
    TestAnalysis & operator=(const TestAnalysis& x);

    //@{
    /// Histograms
    AIDA::IHistogram1D* histTot_;
    AIDA::IHistogram1D* histChTot_;
    AIDA::IHistogram1D* histUnchTot_;
    AIDA::IHistogram1D* histHadrTot_;
    AIDA::IHistogram1D* histHadrChTot_;
    AIDA::IHistogram1D* histHadrUnchTot_;
    //@}

  };

}

#endif
