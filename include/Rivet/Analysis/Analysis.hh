// -*- C++ -*-
#ifndef RIVET_Analysis_HH
#define RIVET_Analysis_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis/Analysis.fhh"
#include "Rivet/Projections/Projection.fhh"
#include "Rivet/Constraints.hh"
#include "Rivet/AnalysisHandler.fhh"
#include "Rivet/Event.fhh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/RivetAIDA.fhh"


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
    
    /// The AnalysisHandler is a friend.
    friend class AnalysisHandler;

  public:

    /// Factory method for getting Analyses.
    static Analysis& getAnalysis(const AnalysisName atype = ANALYSIS_TEST);

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    inline Analysis() 
      : _theHandler(0), _madeHistoDir(false)
    { 
      setBeams(ANY, ANY);
    }

    /// The destructor.
    inline virtual ~Analysis() { }
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

    /// Return the Cuts object of this analysis object. Derived
    /// classes should re-implement this function to return the combined
    /// RivetInfo object of this object and of any Projection objects
    /// upon which this depends.
    virtual const Cuts getCuts() const;

    /// Return the pair of incoming beams required by this analysis.
    inline virtual const BeamPair& getBeams() const {
      return _beams;
    }

    /// Is this analysis compatible with the named quantity set at the supplied value?
    inline virtual const bool isCompatible(const string& quantity, const double value) const {
      Cuts::const_iterator cut = getCuts().find(quantity);
      if (cut == getCuts().end()) return true;
      if (value > cut->second.lowerthan()) return false;
      if (value < cut->second.higherthan()) return false;
      return true;
    }

    /// Is this analysis able to run on the supplied pair of beams?
    inline virtual const bool isCompatible(const ParticleName& beam1, const ParticleName& beam2) const {
      BeamPair beams(beam1, beam2);
      return compatible(beams, getBeams());
      /// @todo Need to also check internal consistency of the analysis' 
      /// beam requirements with those of the projections it uses.
    }

    /// Is this analysis able to run on the BeamPair @a beams ?
    inline virtual const bool isCompatible(const BeamPair& beams) const {
      return compatible(beams, getBeams());
      /// @todo Need to also check internal consistency of the analysis' 
      /// beam requirements with those of the projections it uses.
    }

    /// Access the controlling AnalysisHandler object.
    inline AnalysisHandler& getHandler() const {
      return *_theHandler;
    }

    /// Get the name of the analysis
    inline virtual string getName() const {
      return "";
    }


  protected:
    /// Get a Log object based on the getName() property of the calling analysis object.
    Log& getLog();

    /// Is this analysis able to run on the BeamPair @a beams ?
    virtual const bool checkConsistency() const;
    
    /// Get all the projections used by this analysis, including recursion. 
    /// WARNING: No caching or loop-avoidance is implemented at the moment.
    set<Projection*> getProjections() const;

  protected:
    /// @name AIDA analysis infrastructure.
    //@{
    /// Access the AIDA analysis factory of the controlling AnalysisHandler object.
    AIDA::IAnalysisFactory& analysisFactory();

    /// Access the AIDA tree of the controlling AnalysisHandler object.
    AIDA::ITree& tree();

    /// Access the AIDA histogram factory of the controlling AnalysisHandler object.
    AIDA::IHistogramFactory& histogramFactory();
    //@}

    /// @name Internal histogram booking (for use by Analysis sub-classes).
    //@{

    /// Book a 1D histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const string& name, const string& title, 
                                        const int nbins, const double lower, const double upper);

    /// Book a 1D histogram with non-uniform bins defined by the vector of bin edges @binedges .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const string& name, const string& title, 
                                        const vector<double>& binedges);

    /// Book a 1D histogram based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' getName() property.
    /// @todo Implement auto-binning histo booking methods!
    AIDA::IHistogram1D* bookHistogram1D(const unsigned int datasetId, const unsigned int xAxisId, 
                                        const unsigned int yAxisId, const string& title);


    /// @todo Add profile histograms.


    /// Make the histogram directory.
    void makeHistoDir();


    /// Get the canonical AIDA histogram path for this analysis.
    const string getHistoDir() const {
        return "/" + getName();
    }
    //@}

  protected:

    /// Set the colliding beam pair.
    inline virtual void setBeams(const ParticleName& beam1, const ParticleName& beam2) {
      _beams.first = beam1;
      _beams.second = beam2;
    }

    /// Collection of pointers to projections, for automatically combining constraints.
    set<Projection*> _projections;

  private:

    /// Parameter constraints.
    Cuts _cuts;

    /// Allowed beam-type pair.
    BeamPair _beams;

    /// The controlling AnalysisHandler object.
    AnalysisHandler* _theHandler;

    /// Flag to indicate whether the histogram directory is present
    bool _madeHistoDir;

  private:
    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    Analysis& operator=(const Analysis&);
  };


}

#endif
