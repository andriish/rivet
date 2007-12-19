// -*- C++ -*-
#ifndef RIVET_Analysis_HH
#define RIVET_Analysis_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.fhh"
#include "Rivet/Projection.hh"
#include "Rivet/Constraints.hh"
#include "Rivet/AnalysisHandler.fhh"
#include "Rivet/Event.fhh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// @brief This is the base class of all analysis classes in Rivet.
  ///
  /// There are
  /// three virtual functions which should be implemented in base classes:
  ///
  /// void init() is called by Rivet before a run is started. Here the
  /// analysis class should book necessary histograms. The needed
  /// projections should probably rather be constructed in the
  /// constructor.
  ///
  /// void analyze(const Event&) is called once for each event. Here the
  /// analysis class should apply the necessary Projections and fill the
  /// histograms.
  ///
  /// void finalize() is called after a run is finished. Here the analysis
  /// class should do whatever manipulations are necessary on the
  /// histograms. Writing the histograms to a file is, however, done by
  /// the Rivet class.
  class Analysis {
    
    /// The AnalysisHandler is a friend.
    friend class AnalysisHandler;

  public:

    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor.
    Analysis() 
      : _theHandler(0), _madeHistoDir(false)
    { 
      _gotCrossSection = false;
      setBeams(ANY, ANY);
      setNeedsCrossSection(false);
    }

    /// The destructor.
    virtual ~Analysis() { }
    //@}

  public:

    /// @name Main analysis methods
    //@{
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
    //@}

  public:

    /// @name Metadata
    //@{
    /// Get the name of the analysis. By default this is computed by
    /// combining the results of the experiment, year and Spires ID 
    /// metadata methods and you should only override it if there's a 
    /// good reason why those won't work.
    virtual string getName() const { 
      return getExpt() + "_" + getYear() + "_S" + getSpiresId();
    }

    /// Get a description of the analysis.
    virtual string getSpiresId() const = 0;

    /// Get a description of the analysis.
    virtual string getDescription() const {
      return "";
    }
    
    /// Experiment which performed and published this analysis.
    virtual string getExpt() const = 0;

    /// When published (preprint year according to SPIRES).
    virtual string getYear() const = 0;

    /// Journal, and preprint references.
    virtual vector<string> getReferences() const {
      vector<string> ret;
      return ret;
    }
    //@}

  public:
    /// Return the Cuts object of this analysis object. Derived
    /// classes should re-implement this function to return the combined
    /// RivetInfo object of this object and of any Projection objects
    /// upon which this depends.
    virtual const Cuts getCuts() const;

    /// Return the pair of incoming beams required by this analysis.
    virtual const BeamPair& getBeams() const {
      return _beams;
    }

    /// Is this analysis compatible with the named quantity set at the supplied value?
    virtual const bool isCompatible(const string& quantity, const double value) const {
      Cuts::const_iterator cut = getCuts().find(quantity);
      if (cut == getCuts().end()) return true;
      if (value < cut->second.getHigherThan()) return false;
      if (value > cut->second.getLowerThan()) return false;
      return true;
    }

    /// Is this analysis able to run on the supplied pair of beams?
    virtual const bool isCompatible(const ParticleName& beam1, const ParticleName& beam2) const {
      BeamPair beams(beam1, beam2);
      return compatible(beams, getBeams());
      /// @todo Need to also check internal consistency of the analysis' 
      /// beam requirements with those of the projections it uses.
    }

    /// Is this analysis able to run on the BeamPair @a beams ?
    virtual const bool isCompatible(const BeamPair& beams) const {
      return compatible(beams, getBeams());
      /// @todo Need to also check internal consistency of the analysis' 
      /// beam requirements with those of the projections it uses.
    }

    /// Access the controlling AnalysisHandler object.
    AnalysisHandler& getHandler() const {
      return *_theHandler;
    }

    /// Normalize the given histogram. After this call the histogram
    /// will have been transformed to a DataPointSet with the same
    /// name and path, and the old histogram will be deleted, and the
    /// pointer will be set to zero.
    void normalize(AIDA::IHistogram1D*& histo, const double norm=1.0);

    /// Set the cross section from the generator
    Analysis& setCrossSection(const double& xs) {
      _crossSection = xs;
      _gotCrossSection = true;
      return *this;
    }

    /// Return true if this analysis needs to know the process cross-section.
    bool needsCrossSection() const {
      return _needsCrossSection;
    }
    
  protected:

    /// Get the process cross-section. Throws if this hasn't been set.
    const double& crossSection() {
      if (!_gotCrossSection) {
        string errMsg = "You did not set the cross section for the analysis " + getName();
        throw runtime_error(errMsg);
      }
      return _crossSection;
    }
    
    /// Get a Log object based on the getName() property of the calling analysis object.
    Log& getLog();

    /// Get the number of events seen (via the analysis handler). Use in the
    /// finalize phase only.
    size_t numEvents();

    /// Get the sum of event weights seen (via the analysis handler). Use in the
    /// finalize phase only.
    double sumOfWeights();

    /// Is this analysis able to run on the BeamPair @a beams ?
    virtual const bool checkConsistency() const;
    
    /// Get all the projections used by this analysis, including recursion. 
    /// WARNING: No caching or loop-avoidance is implemented at the moment.
    set<ProjectionPtr> getProjections() const;

  protected:
    /// @name AIDA analysis infrastructure.
    //@{
    /// Access the AIDA analysis factory of the controlling AnalysisHandler object.
    AIDA::IAnalysisFactory& analysisFactory();

    /// Access the AIDA tree of the controlling AnalysisHandler object.
    AIDA::ITree& tree();

    /// Access the AIDA histogram factory of the controlling AnalysisHandler object.
    AIDA::IHistogramFactory& histogramFactory();

    /// Access the AIDA histogram factory of the controlling AnalysisHandler object.
    AIDA::IDataPointSetFactory& datapointsetFactory();

    /// Get the canonical AIDA histogram path for this analysis.
    const string getHistoDir() const {
        return "/" + getName();
    }
    //@}


    /// @name Internal histogram booking (for use by Analysis sub-classes).
    //@{

    /// Book a 1D histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const string& name, const string& title, 
                                        const size_t nbins, const double lower, const double upper);

    /// Book a 1D histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const string& name, const string& title, 
                                        const vector<double>& binedges);

    /// Book a 1D histogram based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' getName() property.
    AIDA::IHistogram1D* bookHistogram1D(const size_t datasetId, const size_t xAxisId, 
                                        const size_t yAxisId, const string& title);
    //@}


    /// @name Internal profile histogram booking (for use by Analysis sub-classes).
    //@{

    /// Book a 1D profile histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IProfile1D* bookProfile1D(const string& name, const string& title, 
                                        const size_t nbins, const double lower, const double upper);

    /// Book a 1D profile histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IProfile1D* bookProfile1D(const string& name, const string& title, 
				        const vector<double>& binedges);

    /// Book a 1D profile histogram based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' getName() property.
    AIDA::IProfile1D* bookProfile1D(const size_t datasetId, const size_t xAxisId, 
                                        const size_t yAxisId, const string& title);
    //@}


    /// @name Internal data point set booking (for use by Analysis sub-classes).
    //@{

    /// Book a 2-dimensional data point set.
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IDataPointSet* bookDataPointSet(const string& name, const string& title);


    /// Book a 2-dimensional data point set with equally spaced points in a range.
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IDataPointSet* bookDataPointSet(const string& name, const string& title, 
                                          const size_t npts, const double lower, const double upper);

    /// Book a 2-dimensional data point set based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings (x-errors) will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' getName() property.
    AIDA::IDataPointSet* bookDataPointSet(const size_t datasetId, const size_t xAxisId, 
                                          const size_t yAxisId, const string& title);
    //@}


  private:

    /// Make the histogram directory.
    void makeHistoDir();


  protected:

    /// Set the colliding beam pair.
    Analysis& setBeams(const ParticleName& beam1, const ParticleName& beam2) {
      _beams.first = beam1;
      _beams.second = beam2;
      return *this;
    }

    /// Add a cut
    Analysis& addCut(const string& quantity, const Comparison& comparison, const double value) {
      _cuts.addCut(quantity, comparison, value);
      return *this;
    }

    /// Add a projection dependency to the projection list.
    Analysis& addProjection(Projection& proj) {
      ProjectionPtr pp(&proj);
      _projections.insert(pp);
      getLog() << Log::DEBUG << " Inserting projection at: " << &proj << endl;
      getLog() << Log::DEBUG << " Inserter/insertee: " << this->getName() << " inserts " << proj.getName() << endl;
      return *this;
    }

    Analysis & setNeedsCrossSection(bool needed){
      _needsCrossSection = needed;
      return *this;
    }
    
    /// Collection of pointers to projections, for automatically combining constraints.
    set<ProjectionPtr> _projections;

  private:

    /// @name Cross-section variables
    //@{
    double _crossSection;
    bool _gotCrossSection;
    bool _needsCrossSection;
    //@}
    
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
