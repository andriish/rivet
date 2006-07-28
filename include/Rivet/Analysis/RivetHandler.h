// -*- C++ -*-
#ifndef RIVET_RivetHandler_H
#define RIVET_RivetHandler_H
//
// This is the declaration of the RivetHandler class.
//

#include "Rivet/Config/Rivet.h"
#include "Rivet/Analysis/AnalysisBase.h"
#include "Rivet/Projections/Event.h"
#include "RivetHandler.fh"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"

namespace Rivet {

/**
 * RivetHandler is a concrete class which handles a number of analysis
 * objects to be applied to generated events.
 */
class RivetHandler {

public:

  /** Typedef a vector of analysis objects. */
  typedef vector<AnaPtr> AnaVector;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The standard constructor.
   * @param filename the name of the file where histograms are stored.
   * @param storetype a string indicating to the AIDA analysis factory
   * how to store the histograms. Which strings are allowed depends on
   * actual AIDA implementation used. To output in standard AIDA XML the
   * string is typically "xml".
   * @param afac an AIDA analysis factory object. The caller must make
   * sure that the lifetime of the factory object exceeds the RivetHandler
   * object.
   */
  RivetHandler(string filename, string storetype,
	       AIDA::IAnalysisFactory & afac);

  /**
   * Default constructor: no AIDA analysis factory required
   */
  RivetHandler();

  /**
   * The destructor is not virtual as this class should not be
   * inherited from.
   */
  ~RivetHandler();
  //@}

public:

  /**
   * Add an object of base class AnalysisBase to the list of analysis
   * objects to be used in a run. Note that the object will be copied,
   * and the agrument object given will not be used and can be
   * discarded directly.
   */
  template <typename ANA>
  inline void addAnalysis(const ANA & ana);

  /**
   * Initialize a run. If this run is to be joined together with other
   * runs, \a N should be set to the total number of runs to be
   * combined, and \a i should be the index of this run. This function
   * will initialize the histogram factory and then call the
   * AnalysisBase::init() function of all included analysis objects.
   */
  void init(int i = 0, int N = 0);

  /**
   * Analyze the given \a event. This function will call the
   * AnalysisBase::analyze() function of all included analysis objects.
   */
  void analyze(const GenEvent & event);

  /**
   * Finalize a run. This function first calls the
   * AnalysisBase::finalize() functions of all included analysis
   * objects and then writes out all histograms to a file.
   */
  void finalize();

  /**
   * Return a RivetInfo object containing all parameters of all
   * included analysis handlers and all their included Projections.
   */
  RivetInfo info() const;

  /**
   * The AIDA analysis factory.
   */
  inline AIDA::IAnalysisFactory & analysisFactory();

  /**
   * The AIDA tree object.
   */
  inline AIDA::ITree & tree();

  /**
   * The AIDA histogram factory.
   */
  inline AIDA::IHistogramFactory & histogramFactory();


private:

  /**
   * The vector of AnalysisBase objects to be used.
   */
  AnaVector anaVector;

  /**
   * If non-zero the number of runs to be combined into one analysis.
   */
  int nRun;

  /**
   * If non-zero, the index of this run, if a part of several runs to
   * be combined into one.
   */
  int iRun;

  /**
   * The AIDA analysis factory.
   */
  AIDA::IAnalysisFactory * theAnalysisFactory;

  /**
   * The AIDA tree object.
   */
  AIDA::ITree * theTree;

  /**
   * The AIDA histogram factory.
   */
  AIDA::IHistogramFactory * theHistogramFactory;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RivetHandler & operator=(const RivetHandler &);

  /**
   * The copy constructor is private and must never be called.  In
   * fact, it should not even be implemented.
   */
  inline RivetHandler(const RivetHandler &);

  /*
   * *** ATTENTION *** Here we should have a histogram factory.
   */

};

}

#include "Rivet/Analysis/RivetHandler.icc"

#endif /* RIVET_RivetHandler_H */
