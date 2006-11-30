// -*- C++ -*-
#ifndef RIVET_TestAnalysis_H
#define RIVET_TestAnalysis_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"

namespace Rivet {

  /// This class just measures the charged multiplicity
  class TestAnalysis : public Analysis {

  public:

    /// Default constructor.
    inline TestAnalysis();

    /// Copy constructor.
    inline TestAnalysis(const TestAnalysis &);

    /// Destructor
    virtual ~TestAnalysis();

  public:

    virtual void init();

    virtual void analyze(const Event & event);

    virtual void finalize();

    /// Return the RivetInfo object of this analysis object.
    virtual RivetInfo getInfo() const;

  private:

    /// The Multiplicity projector used by this analysis.
    Multiplicity mult;

  private:

    TestAnalysis & operator=(const TestAnalysis &);

  };

}

#include "Rivet/Analysis/TestAnalysis.icc"

#endif
