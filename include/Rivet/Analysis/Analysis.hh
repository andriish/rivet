// -*- C++ -*-
#ifndef RIVET_Analysis_HH
#define RIVET_Analysis_HH
//
// This is the declaration of the Analysis base class.
//

#include "Rivet/Analysis/Analysis.fhh"

#include "Rivet/RivetInfo.hh"
#include "Rivet/RivetHandler.fhh"
#include "Rivet/Tools/Event/Event.fhh"


namespace AIDA {
  class IAnalysisFactory;
  class IHistogramFactory;
  class ITree;
  class IHistogram1D;
}


namespace Rivet {

  /// This is the base class of all analysis classes in Rivet. There are
  /// three virtual functions which should be implemented in base classes:
  ///
  /// void init() is called by Rivet before a run is started. Here the
  /// analysis class should book necessary histograms. The needed
  /// projections should probably rather be constructed in the
  /// constructor.
  ///
  /// void analyze(const Event &) is called once for each event. Here the
  /// analysis class should apply the necessary Projections and fill the
  /// histograms.
  ///
  /// void finish() is called after a run is finished. Here the analysis
  /// class should do whatever manipulations are necessary on the
  /// histograms. Writing the histograms to a file is, however, done by
  /// the Rivet class.
  class Analysis {
    
    /// The RivetHandler is a friend.
    friend class RivetHandler;

  public:

    /// Factory method for getting Analyses.
    static Analysis& getAnalysis(const AnalysisName atype = ANALYSIS_TEST);

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    inline Analysis() 
      : theHandler(0), _madeHistoDir(false)
    { }

    /// The copy constructor.
    inline Analysis(const Analysis& x) 
      : theHandler(x.theHandler), info(x.info), _madeHistoDir(false) 
    { }

    /// The destructor.
    virtual ~Analysis();
    //@}

  public:

    /// Initialize this analysis object. A concrete class should here
    /// book all necessary histograms. An overridden function must make
    /// sure it first calls the base class function.
    virtual void init() = 0;

    /// Analyze one event. A concrete class should here apply the
    /// necessary projections on the \a event and fill the relevant
    /// histograms. An overridden function must make sure it first calls
    /// the base class function.
    virtual void analyze(const Event& event) = 0;

    /// Finalize this analysis object. A concrete class should here make
    /// all necessary operations on the histograms. Writing the
    /// histograms to a file is, however, done by the Rivet class. An
    /// overridden function must make sure it first calls the base class
    /// function.
    virtual void finalize() = 0;

    /// Return the RivetInfo object of this analysis object. Derived
    /// classes should re-implement this function to return the combined
    /// RivetInfo object of this object and of any Projection objects
    /// upon which this depends.
    virtual RivetInfo getInfo() const;

    /// Access the controlling RivetHandler object.
    inline RivetHandler& handler() const {
      return *theHandler;
    }

    /// Get the name of the analysis
    virtual std::string name() const {
      return "";
    }


  protected:

    /// Access the AIDA analysis factory of the controlling RivetHandler object.
    AIDA::IAnalysisFactory& analysisFactory();

    /// Access the AIDA tree of the controlling RivetHandler object.
    AIDA::ITree& tree();

    /// Access the AIDA histogram factory of the controlling RivetHandler object.
    AIDA::IHistogramFactory& histogramFactory();

    /// Book a 1D histogram (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const std::string& name, const std::string& title, 
                                        const int nbins, const double lower, const double upper);

    /// Make the histogram directory
    void makeHistoDir();

    const std::string histoDir() const {
        return "/" + name();
    }

  private:

    /// The controlling RivetHandler object.
    RivetHandler* theHandler;

    /// The object containing the parameters of this analysis object to
    /// be communicated to the outside world.
    RivetInfo info;

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    Analysis& operator=(const Analysis&);

    /// Flag to indicate whether the histogram directory is present
    bool _madeHistoDir;

  };


}

#endif /* RIVET_Analysis_H */
