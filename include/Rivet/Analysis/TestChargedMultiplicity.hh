// -*- C++ -*-
#ifndef RIVET_TestChargedMultiplicity_H
#define RIVET_TestChargedMultiplicity_H

#include "Rivet/Analysis/AnalysisBase.hh"
#include "Rivet/Projections/Multiplicity.hh"

namespace Rivet {

/// This class just measures the charged multiplicity
class TestChargedMultiplicity : public AnalysisBase {

public:

  /// Default constructor.
  inline TestChargedMultiplicity();

  /// Copy constructor.
  inline TestChargedMultiplicity(const TestChargedMultiplicity &);

  /// Destructor
  virtual ~TestChargedMultiplicity();

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

  TestChargedMultiplicity & operator=(const TestChargedMultiplicity &);

};


inline TestChargedMultiplicity::TestChargedMultiplicity()
  : mult() {}

inline TestChargedMultiplicity::TestChargedMultiplicity(const TestChargedMultiplicity & x)
  : AnalysisBase(x), mult(x.mult) {}


}

#endif
