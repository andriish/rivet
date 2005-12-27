// -*- C++ -*-
#ifndef RIVET_TestMultiplicity_H
#define RIVET_TestMultiplicity_H
//
// This is the declaration of the TestMultiplicity class.
//

#include "Rivet/Analysis/AnalysisBase.h"
#include "Rivet/Projections/FinalStateHCM.h"
#include "Rivet/Projections/FinalStateProjection.h"

namespace Rivet {

/**
 * This class simply measures the total multiplicity. It is only
 * intended for testing purposes.
 */
class TestMultiplicity: public AnalysisBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline TestMultiplicity();

  /**
   * The copy constructor.
   */
  inline TestMultiplicity(const TestMultiplicity &);

  /**
   * The destructor.
   */
  virtual ~TestMultiplicity();
  //@}

public:

  /**
   * Initialize this analysis object. A concrete class should here
   * book all necessary histograms. An overridden function must make
   * sure it first calls the base class function.
   */
  virtual void init();

  /**
   * Analyze one event. A concrete class should here apply the
   * necessary projections on the \a event and fill the relevant
   * histograms. An overridden function must make sure it first calls
   * the base class function.
   */
  virtual void analyze(const Event & event);

  /**
   * Finalize this analysis object. A concrete class should here make
   * all necessary operations on the histograms. Writing the
   * histograms to a file is, however, done by the Rivet class. An
   * overridden function must make sure it first calls the base class
   * function.
   */
  virtual void finalize();

  /**
   * Return the RivetInfo object of this analysis object. Derived
   * classes should re-implement this function to return the combined
   * RivetInfo object of this object and of any Projection objects
   * upon which this depends.
   */
  virtual RivetInfo getInfo() const;

private:

  /**
   * The FinalStateProjector used.
   */
  FinalStateProjection fsproj;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TestMultiplicity & operator=(const TestMultiplicity &);

};

}

#include "TestMultiplicity.icc"

#endif /* RIVET_TestMultiplicity_H */
