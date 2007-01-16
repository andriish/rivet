// -*- C++ -*-
#ifndef RIVET_Analysis_H
#define RIVET_Analysis_H
//
// This is the declaration of the Analysis base class.
//

#include "Rivet/RivetInfo.hh"
#include "Rivet/RivetHandler.fhh"
#include "Rivet/Tools/Event/Event.fhh"
#include "Analysis.fhh"

namespace AIDA {
  class IAnalysisFactory;
  class IHistogramFactory;
  class ITree;
}


namespace Rivet {

  /**
   * This is the base class of all analysis classes in Rivet. There are
   * three virtual functions which should be implemented in base classes:
   *
   * void init() is called by Rivet before a run is started. Here the
   * analysis class should book necessary histograms. The needed
   * projections should probably rather be constructed in the
   * constructor.
   *
   * void analyze(const Event &) is called once for each event. Here the
   * analysis class should apply the necessary Projections and fill the
   * histograms.
   *
   * void finish() is called after a run is finsished. Here the analysis
   * class should do whatever manipulations are necessay on the
   * histograms. Writing the histograms to a file is, however, done by
   * the Rivet class.
   */
  class Analysis {
    
    /**
     * The RivetHandler is a friend;
     */
    friend class RivetHandler;

  public:

    /// Factory method for getting Analyses
    static Analysis& getAnalysis(const AnalysisName atype = ANALYSIS_TEST);

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * The default constructor.
     */
    inline Analysis();

    /**
     * The copy constructor.
     */
    inline Analysis(const Analysis &);

    /**
     * The destructor.
     */
    virtual ~Analysis();
    //@}

  public:

    /**
     * Initialize this analysis object. A concrete class should here
     * book all necessary histograms. An overridden function must make
     * sure it first calls the base class function.
     */
    virtual void init() = 0;

    /**
     * Analyze one event. A concrete class should here apply the
     * necessary projections on the \a event and fill the relevant
     * histograms. An overridden function must make sure it first calls
     * the base class function.
     */
    virtual void analyze(const Event & event) = 0;

    /**
     * Finalize this analysis object. A concrete class should here make
     * all necessary operations on the histograms. Writing the
     * histograms to a file is, however, done by the Rivet class. An
     * overridden function must make sure it first calls the base class
     * function.
     */
    virtual void finalize() = 0;

    /**
     * Return the RivetInfo object of this analysis object. Derived
     * classes should re-implement this function to return the combined
     * RivetInfo object of this object and of any Projection objects
     * upon which this depends.
     */
    virtual RivetInfo getInfo() const;

    /**
     * Access the controlling RivetHandler object.
     */
    inline RivetHandler & handler();

    /**
     * Access the AIDA analysis factory of the controlling RivetHandler
     * object.
     */
    AIDA::IAnalysisFactory & analysisFactory();

    /**
     * Access the AIDA tree of the controlling RivetHandler object.
     */
    AIDA::ITree & tree();

    /**
     * Access the AIDA histogram factory of the controlling RivetHandler
     * object.
     */
    AIDA::IHistogramFactory & histogramFactory();

  private:

    /**
     * The controlling RivetHandler object.
     */
    RivetHandler * theHandler;

    /**
     * The object containing the parameters of this analysis object to
     * be communicated to the outside world.
     */
    RivetInfo info;

  private:

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    Analysis & operator=(const Analysis &);

  };

}

#include "Rivet/Analysis/Analysis.icc"

#endif /* RIVET_Analysis_H */
