// -*- C++ -*-
#ifndef RIVET_Analysis_HH
#define RIVET_Analysis_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.fhh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/Constraints.hh"
#include "Rivet/AnalysisHandler.fhh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/RivetAIDA.fhh"

/// \def vetoEvent
/// Preprocessor define for vetoing events, including the log message and return.
#define vetoEvent(E) { vetoEventWeight(E); getLog() << Log::DEBUG << "Vetoing event on line " << __LINE__ << " of " << __FILE__ << endl; return; }


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
  class Analysis : public ProjectionApplier {
    
    /// The AnalysisHandler is a friend.
    friend class AnalysisHandler;

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    Analysis();

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
    /// Metadata is used for querying from the command line and also for
    /// building web pages and the analysis pages in the Rivet manual.
    //@{
    /// Get the name of the analysis. By default this is computed by
    /// combining the results of the experiment, year and Spires ID 
    /// metadata methods and you should only override it if there's a 
    /// good reason why those won't work.
    virtual std::string name() const {
      return experiment() + "_" + year() + "_S" + spiresId();
    }

    /// Get a description of the analysis.
    virtual std::string spiresId() const = 0;

    /// @brief Names & emails of paper/analysis authors.
    /// Names and email of authors in 'NAME <EMAIL>' format. The first
    /// name in the list should be the primary contact person.
    //MD virtual std::vector<std::string> authors() const = 0;

    /// @brief Get a short description of the analysis.
    /// Short (one sentence) description used as an index entry.
    /// Use @a description() to provide full descriptive paragraphs
    /// of analysis details.
    virtual std::string summary() const = 0;

    /// @brief Get a full description of the analysis.
    /// Full textual description of this analysis, what it is useful for,
    /// what experimental techniques are applied, etc. Should be treated
    /// as a chunk of restructuredText (http://docutils.sourceforge.net/rst.html),
    /// with equations to be rendered as LaTeX with amsmath operators.
    //MD virtual std::string description() const = 0;

    /// @brief Information about the events needed as input for this analysis.
    /// Event types, energies, kinematic cuts, particles to be considered 
    /// stable, etc. etc. Should be treated as a restructuredText bullet list
    /// (http://docutils.sourceforge.net/rst.html)
    //MD virtual std::string runInfo() const = 0;
    
    /// Experiment which performed and published this analysis.
    virtual std::string experiment() const = 0;

    /// Collider on which the experiment ran.
    //MD virtual string collider() const = 0;

    /// Incoming beams required by this analysis.
    virtual const BeamPair& beams() const {
      return _beams;
    }

    /// @brief When the original experimental analysis was published.
    /// When the refereed paper on which this is based was published, 
    /// according to SPIRES.
    virtual std::string year() const = 0;

    /// Journal, and preprint references.
    //MD virtual std::vector<std::string> references() const = 0;

    /// Whether this analysis is trusted (in any way!)
    virtual std::string status() const {
      return "UNVALIDATED";
    }
    //@}


  public:

    /// Return the pair of incoming beams required by this analysis.
    virtual const BeamPair& requiredBeams() const {
      return _beams;
    }

    /// Is this analysis able to run on the supplied pair of beams?
    virtual const bool isCompatible(const ParticleName& beam1, const ParticleName& beam2) const {
      BeamPair beams(beam1, beam2);
      return compatible(beams, requiredBeams());
      /// @todo Need to also check internal consistency of the analysis' 
      /// beam requirements with those of the projections it uses.
    }

    /// Is this analysis able to run on the BeamPair @a beams ?
    virtual const bool isCompatible(const BeamPair& beams) const {
      return compatible(beams, requiredBeams());
      /// @todo Need to also check internal consistency of the analysis' 
      /// beam requirements with those of the projections it uses.
    }

    /// Access the controlling AnalysisHandler object.
    AnalysisHandler& handler() const {
      return *_analysishandler;
    }

    /// Normalize the given histogram, @a histo. After this call the 
    /// histogram will have been transformed to a DataPointSet with the 
    /// same name and path. It has the same effect as 
    /// @c scale(histo, norm/sumOfWeights).
    /// @param histo The histogram to be normalised.
    /// @param norm The new area of the histogram.
    /// @warning The old histogram will be deleted, and its pointer set to zero.
    void normalize(AIDA::IHistogram1D*& histo, const double norm=1.0);

    /// Multiplicatively scale the given histogram, @a histo. After this call the 
    /// histogram will have been transformed to a DataPointSet with the same name and path.  
    /// @param histo The histogram to be scaled.
    /// @param scale The factor used to multiply the histogram bin heights.
    /// @warning The old histogram will be deleted, and its pointer set to zero.
    void scale(AIDA::IHistogram1D*& histo, const double scale);

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

    /// Get the process cross-section in pb. Throws if this hasn't been set.
    const double& crossSection() const {
      if (!_gotCrossSection) {
        string errMsg = "You did not set the cross section for the analysis " + name();
        throw Error(errMsg);
      }
      return _crossSection;
    }
    
    // Veto the current event by subtracting its weight from the tally.
    void vetoEventWeight(const Event& e) {
      _vetoedWeightSum += e.weight();
    }

    /// Get a Log object based on the name() property of the calling analysis object.
    Log& getLog() const;

    /// Get the number of events seen (via the analysis handler). Use in the
    /// finalize phase only.
    size_t numEvents() const;

    /// Get the sum of event weights seen (via the analysis handler). Use in the
    /// finalize phase only.
    double sumOfWeights() const;
    

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
    const std::string histoDir() const {
        return "/" + name();
    }
    //@}


    /// @name Internal histogram booking (for use by Analysis sub-classes).
    //@{

    /// Book a 1D histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const std::string& name, const std::string& title, 
                                        const size_t nbins, const double lower, const double upper);

    /// Book a 1D histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const std::string& name, const std::string& title, 
                                        const std::vector<double>& binedges);

    /// Book a 1D histogram based on the name in the corresponding AIDA
    /// file. The binnings will be obtained by reading the bundled AIDA data
    /// record file with the same filename as the analysis' name() property.
    AIDA::IHistogram1D* bookHistogram1D(const std::string& name, const std::string& title);

    /// Book a 1D histogram based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' name() property.
    AIDA::IHistogram1D* bookHistogram1D(const size_t datasetId, const size_t xAxisId, 
                                        const size_t yAxisId, const std::string& title);
    //@}


    /// @name Internal profile histogram booking (for use by Analysis sub-classes).
    //@{

    /// Book a 1D profile histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IProfile1D* bookProfile1D(const std::string& name, const std::string& title, 
                                        const size_t nbins, const double lower, const double upper);

    /// Book a 1D profile histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IProfile1D* bookProfile1D(const std::string& name, const std::string& title, 
                                    const std::vector<double>& binedges);

    /// Book a 1D profile histogram based on the name in the corresponding AIDA
    /// file. The binnings will be obtained by reading the bundled AIDA data
    /// record file with the same filename as the analysis' name() property.
    AIDA::IProfile1D* bookProfile1D(const std::string& name, const std::string& title);

    /// Book a 1D profile histogram based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' name() property.
    AIDA::IProfile1D* bookProfile1D(const size_t datasetId, const size_t xAxisId, 
                                        const size_t yAxisId, const std::string& title);
    //@}


    /// @name Internal data point set booking (for use by Analysis sub-classes).
    //@{

    /// Book a 2-dimensional data point set.
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IDataPointSet* bookDataPointSet(const std::string& name, const std::string& title);


    /// Book a 2-dimensional data point set with equally spaced points in a range.
    /// (NB. this returns a pointer rather than a reference since it will 
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly 
    /// get the pointer from a reference before they can use it!)
    AIDA::IDataPointSet* bookDataPointSet(const std::string& name, const std::string& title, 
                                          const size_t npts, const double lower, const double upper);

    /// Book a 2-dimensional data point set based on the corresponding AIDA data
    /// file. The binnings (x-errors) will be obtained by reading the bundled
    /// AIDA data record file of the same filename as the analysis' name()
    /// property.
    //AIDA::IDataPointSet* bookDataPointSet(const std::string& name, const std::string& title);

    /// Book a 2-dimensional data point set based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings (x-errors) will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' name() property.
    AIDA::IDataPointSet* bookDataPointSet(const size_t datasetId, const size_t xAxisId, 
                                          const size_t yAxisId, const std::string& title);
    //@}


  private:

    /// @name Utility functions
    //@{

    /// Make the histogram directory.
    void _makeHistoDir();

    /// Get the bin edges for this paper from the reference AIDA file, and cache them.
    void _cacheBinEdges();

    /// Get the x-axis points for this paper from the reference AIDA file, and cache them.
    void _cacheXAxisData();

    /// Make the axis code string (dsDD-xXX-yYY)
    string _makeAxisCode(const size_t datasetId, const size_t xAxisId, const size_t yAxisId) const;
    //@}


  protected:

    /// Set the colliding beam pair.
    Analysis& setBeams(const ParticleName& beam1, const ParticleName& beam2) {
      _beams.first = beam1;
      _beams.second = beam2;
      return *this;
    }

    /// Declare whether this analysis needs to know the process cross-section from the generator.
    Analysis& setNeedsCrossSection(bool needed) {
      _needsCrossSection = needed;
      return *this;
    }
    
  private:

    /// @name Cross-section variables
    //@{
    double _crossSection;
    bool _gotCrossSection;
    bool _needsCrossSection;
    //@}
    
    /// Allowed beam-type pair.
    BeamPair _beams;

    /// The controlling AnalysisHandler object.
    AnalysisHandler* _analysishandler;

    /// Flag to indicate whether the histogram directory is present
    bool _madeHistoDir;

    double _vetoedWeightSum;

    /// Collection of x-axis point data to speed up many autobookings: the 
    /// AIDA reference file should only be read once.
    map<string, vector<DPSXPoint> > _dpsData;

    /// Collection of cached bin edges to speed up many autobookings: the 
    /// AIDA reference file should only be read once.
    map<string, BinEdges> _histBinEdges;

  private:
    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    Analysis& operator=(const Analysis&);
  };

}

#endif
