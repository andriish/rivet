// -*- C++ -*-
#ifndef RIVET_Analysis_HH
#define RIVET_Analysis_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleUtils.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Tools/RivetMT2.hh"
#include "Rivet/Tools/RivetYODA.hh"


/// @def vetoEvent
/// Preprocessor define for vetoing events, including the log message and return.
#define vetoEvent                                                       \
  do { MSG_DEBUG("Vetoing event on line " << __LINE__ << " of " << __FILE__); return; } while(0)


namespace Rivet {

  // convenience for analysis writers
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::stringstream;
  using std::swap;
  using std::numeric_limits;

  // Forward declaration
  class AnalysisHandler;

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

    // /// The default constructor.
    // Analysis();

    /// Constructor
    Analysis(const std::string& name);

    /// The destructor.
    virtual ~Analysis() {}

    //@}


  public:

    /// @name Main analysis methods
    //@{

    /// Initialize this analysis object. A concrete class should here
    /// book all necessary histograms. An overridden function must make
    /// sure it first calls the base class function.
    virtual void init() { }

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
    virtual void finalize() { }

    //@}


  public:

    /// @name Metadata
    /// Metadata is used for querying from the command line and also for
    /// building web pages and the analysis pages in the Rivet manual.
    //@{

    /// Get the actual AnalysisInfo object in which all this metadata is stored.
    const AnalysisInfo& info() const {
      assert(_info && "No AnalysisInfo object :O");
      return *_info;
    }

    /// @brief Get the name of the analysis.
    ///
    /// By default this is computed by combining the results of the experiment,
    /// year and Spires ID metadata methods and you should only override it if
    /// there's a good reason why those won't work.
    virtual std::string name() const {
      return (info().name().empty()) ? _defaultname : info().name();
    }

    /// Get the Inspire ID code for this analysis.
    virtual std::string inspireId() const {
      return info().inspireId();
    }

    /// Get the SPIRES ID code for this analysis (~deprecated).
    virtual std::string spiresId() const {
      return info().spiresId();
    }

    /// @brief Names & emails of paper/analysis authors.
    ///
    /// Names and email of authors in 'NAME \<EMAIL\>' format. The first
    /// name in the list should be the primary contact person.
    virtual std::vector<std::string> authors() const {
      return info().authors();
    }

    /// @brief Get a short description of the analysis.
    ///
    /// Short (one sentence) description used as an index entry.
    /// Use @a description() to provide full descriptive paragraphs
    /// of analysis details.
    virtual std::string summary() const {
      return info().summary();
    }

    /// @brief Get a full description of the analysis.
    ///
    /// Full textual description of this analysis, what it is useful for,
    /// what experimental techniques are applied, etc. Should be treated
    /// as a chunk of restructuredText (http://docutils.sourceforge.net/rst.html),
    /// with equations to be rendered as LaTeX with amsmath operators.
    virtual std::string description() const {
      return info().description();
    }

    /// @brief Information about the events needed as input for this analysis.
    ///
    /// Event types, energies, kinematic cuts, particles to be considered
    /// stable, etc. etc. Should be treated as a restructuredText bullet list
    /// (http://docutils.sourceforge.net/rst.html)
    virtual std::string runInfo() const {
      return info().runInfo();
    }

    /// Experiment which performed and published this analysis.
    virtual std::string experiment() const {
      return info().experiment();
    }

    /// Collider on which the experiment ran.
    virtual std::string collider() const {
      return info().collider();
    }

    /// When the original experimental analysis was published.
    virtual std::string year() const {
      return info().year();
    }

    /// The luminosity in inverse femtobarn
    virtual std::string luminosityfb() const {
      return info().luminosityfb();
    }

    /// Journal, and preprint references.
    virtual std::vector<std::string> references() const {
      return info().references();
    }

    /// BibTeX citation key for this article.
    virtual std::string bibKey() const {
      return info().bibKey();
    }

    /// BibTeX citation entry for this article.
    virtual std::string bibTeX() const {
      return info().bibTeX();
    }

    /// Whether this analysis is trusted (in any way!)
    virtual std::string status() const {
      return (info().status().empty()) ? "UNVALIDATED" : info().status();
    }

    /// Any work to be done on this analysis.
    virtual std::vector<std::string> todos() const {
      return info().todos();
    }


    /// Return the allowed pairs of incoming beams required by this analysis.
    virtual const std::vector<PdgIdPair>& requiredBeams() const {
      return info().beams();
    }
    /// Declare the allowed pairs of incoming beams required by this analysis.
    virtual Analysis& setRequiredBeams(const std::vector<PdgIdPair>& requiredBeams) {
      info().setBeams(requiredBeams);
      return *this;
    }


    /// Sets of valid beam energy pairs, in GeV
    virtual const std::vector<std::pair<double, double> >& requiredEnergies() const {
      return info().energies();
    }

    /// Get vector of analysis keywords
    virtual const std::vector<std::string> & keywords() const {
      return info().keywords();
    }

    /// Declare the list of valid beam energy pairs, in GeV
    virtual Analysis& setRequiredEnergies(const std::vector<std::pair<double, double> >& requiredEnergies) {
      info().setEnergies(requiredEnergies);
      return *this;
    }


    //@}


    /// @name Internal metadata modifying methods
    //@{

    /// Get the actual AnalysisInfo object in which all this metadata is stored (non-const).
    AnalysisInfo& info() {
      assert(_info && "No AnalysisInfo object :O");
      return *_info;
    }

    //@}


    /// @name Run conditions
    //@{

    /// Incoming beams for this run
    const ParticlePair& beams() const;

    /// Incoming beam IDs for this run
    const PdgIdPair beamIds() const;

    /// Centre of mass energy for this run
    double sqrtS() const;

    //@}


    /// @name Analysis / beam compatibility testing
    //@{

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(const ParticlePair& beams) const;

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(PdgId beam1, PdgId beam2, double e1, double e2) const;

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(const PdgIdPair& beams, const std::pair<double,double>& energies) const;

    //@}

    /// Access the controlling AnalysisHandler object.
    AnalysisHandler& handler() const { return *_analysishandler; }


  protected:

    /// Get a Log object based on the name() property of the calling analysis object.
    Log& getLog() const;

    /// Get the process cross-section in pb. Throws if this hasn't been set.
    double crossSection() const;

    /// Get the process cross-section per generated event in pb. Throws if this
    /// hasn't been set.
    double crossSectionPerEvent() const;

    /// @brief Get the number of events seen (via the analysis handler).
    ///
    /// @note Use in the finalize phase only.
    size_t numEvents() const;

    /// @brief Get the sum of event weights seen (via the analysis handler).
    ///
    /// @note Use in the finalize phase only.
    double sumOfWeights() const;

  protected:

    /// @name Histogram paths
    //@{

    /// Get the canonical histogram "directory" path for this analysis.
    const std::string histoDir() const;

    /// Get the canonical histogram path for the named histogram in this analysis.
    const std::string histoPath(const std::string& hname) const;

    /// Get the canonical histogram path for the numbered histogram in this analysis.
    const std::string histoPath(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const;

    /// Get the internal histogram name for given d, x and y (cf. HepData)
    const std::string mkAxisCode(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const;
    /// Alias
    /// @deprecated Prefer the "mk" form, consistent with other "making function" names
    const std::string makeAxisCode(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
      return mkAxisCode(datasetId, xAxisId, yAxisId);
    }

    //@}


    /// @name Histogram reference data
    //@{

    /// Get reference data for a named histo
    /// @todo SFINAE to ensure that the type inherits from YODA::AnalysisObject?
    template <typename T=YODA::Scatter2D>
    const T& refData(const string& hname) const {
      _cacheRefData();
      MSG_TRACE("Using histo bin edges for " << name() << ":" << hname);
      if (!_refdata[hname]) {
        MSG_ERROR("Can't find reference histogram " << hname);
        throw Exception("Reference data " + hname + " not found.");
      }
      return dynamic_cast<T&>(*_refdata[hname]);
    }

    /// Get reference data for a numbered histo
    /// @todo SFINAE to ensure that the type inherits from YODA::AnalysisObject?
    template <typename T=YODA::Scatter2D>
    const T& refData(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
      const string hname = makeAxisCode(datasetId, xAxisId, yAxisId);
      return refData(hname);
    }

    //@}


    /// @name Counter booking
    //@{

    /// Book a counter.
    CounterPtr & book(CounterPtr &, const std::string& name,
                           const std::string& title="");
                           // const std::string& valtitle=""

    /// Book a counter, using a path generated from the dataset and axis ID codes
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    CounterPtr & book(CounterPtr &, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                           const std::string& title="");
                           // const std::string& valtitle=""

    //@}


    /// @name 1D histogram booking
    //@{

    /// Book a 1D histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    Histo1DPtr & book(Histo1DPtr &,const std::string& name,
                           size_t nbins, double lower, double upper,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="");

    /// Book a 1D histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    Histo1DPtr & book(Histo1DPtr &,const std::string& name,
                           const std::vector<double>& binedges,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="");

    /// Book a 1D histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    Histo1DPtr & book(Histo1DPtr &,const std::string& name,
                           const std::initializer_list<double>& binedges,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="");

    /// Book a 1D histogram with binning from a reference scatter.
    Histo1DPtr & book(Histo1DPtr &,const std::string& name,
                           const Scatter2D& refscatter,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="");

    /// Book a 1D histogram, using the binnings in the reference data histogram.
    Histo1DPtr & book(Histo1DPtr &,const std::string& name,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="");

    /// Book a 1D histogram, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    Histo1DPtr & book(Histo1DPtr &,unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="");
    //@}


    /// @name 2D histogram booking
    //@{

    /// Book a 2D histogram with @a nxbins and @a nybins uniformly
    /// distributed across the ranges @a xlower - @a xupper and @a
    /// ylower - @a yupper respectively along the x- and y-axis.
    Histo2DPtr & book(Histo2DPtr &,const std::string& name,
                           size_t nxbins, double xlower, double xupper,
                           size_t nybins, double ylower, double yupper,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="",
                           const std::string& ztitle="");

    /// Book a 2D histogram with non-uniform bins defined by the
    /// vectors of bin edges @a xbinedges and @a ybinedges.
    Histo2DPtr & book(Histo2DPtr &,const std::string& name,
                           const std::vector<double>& xbinedges,
                           const std::vector<double>& ybinedges,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="",
                           const std::string& ztitle="");

    /// Book a 2D histogram with non-uniform bins defined by the
    /// vectors of bin edges @a xbinedges and @a ybinedges.
    Histo2DPtr & book(Histo2DPtr &,const std::string& name,
                           const std::initializer_list<double>& xbinedges,
                           const std::initializer_list<double>& ybinedges,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="",
                           const std::string& ztitle="");

    /// Book a 2D histogram with binning from a reference scatter.
    Histo2DPtr & book(Histo2DPtr &,const std::string& name,
                           const Scatter3D& refscatter,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="",
                           const std::string& ztitle="");

    /// Book a 2D histogram, using the binnings in the reference data histogram.
    Histo2DPtr & book(Histo2DPtr &,const std::string& name,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="",
                           const std::string& ztitle="");

    /// Book a 2D histogram, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    Histo2DPtr & book(Histo2DPtr &,unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                           const std::string& title="",
                           const std::string& xtitle="",
                           const std::string& ytitle="",
                           const std::string& ztitle="");

    //@}


    /// @name 1D profile histogram booking
    //@{

    /// Book a 1D profile histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    Profile1DPtr & book(Profile1DPtr &,  const std::string& name,
                               size_t nbins, double lower, double upper,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    /// Book a 1D profile histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    Profile1DPtr & book(Profile1DPtr &,  const std::string& name,
                               const std::vector<double>& binedges,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    /// Book a 1D profile histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    Profile1DPtr & book(Profile1DPtr &,  const std::string& name,
                               const std::initializer_list<double>& binedges,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    /// Book a 1D profile histogram with binning from a reference scatter.
    Profile1DPtr & book(Profile1DPtr &,  const std::string& name,
                               const Scatter2D& refscatter,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    /// Book a 1D profile histogram, using the binnings in the reference data histogram.
    Profile1DPtr & book(Profile1DPtr &,  const std::string& name,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    /// Book a 1D profile histogram, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    Profile1DPtr & book(Profile1DPtr &,  unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    //@}


    /// @name 2D profile histogram booking
    //@{

    /// Book a 2D profile histogram with @a nxbins and @a nybins uniformly
    /// distributed across the ranges @a xlower - @a xupper and @a ylower - @a
    /// yupper respectively along the x- and y-axis.
    Profile2DPtr & book(Profile2DPtr &,  const std::string& name,
                               size_t nxbins, double xlower, double xupper,
                               size_t nybins, double ylower, double yupper,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="",
                               const std::string& ztitle="");

    /// Book a 2D profile histogram with non-uniform bins defined by the vectorx
    /// of bin edges @a xbinedges and @a ybinedges.
    Profile2DPtr & book(Profile2DPtr &,  const std::string& name,
                               const std::vector<double>& xbinedges,
                               const std::vector<double>& ybinedges,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="",
                               const std::string& ztitle="");

    /// Book a 2D profile histogram with non-uniform bins defined by the vectorx
    /// of bin edges @a xbinedges and @a ybinedges.
    Profile2DPtr & book(Profile2DPtr &,  const std::string& name,
                               const std::initializer_list<double>& xbinedges,
                               const std::initializer_list<double>& ybinedges,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="",
                               const std::string& ztitle="");

    /// Book a 2D profile histogram with binning from a reference scatter.
    // Profile2DPtr bookProfile2D(const std::string& name,
    //                            const Scatter3D& refscatter,
    //                            const std::string& title="",
    //                            const std::string& xtitle="",
    //                            const std::string& ytitle="",
    //                            const std::string& ztitle="");

    // /// Book a 2D profile histogram, using the binnings in the reference data histogram.
    // Profile2DPtr bookProfile2D(const std::string& name,
    //                            const std::string& title="",
    //                            const std::string& xtitle="",
    //                            const std::string& ytitle="",
    //                            const std::string& ztitle="");

    // /// Book a 2D profile histogram, using the binnings in the reference data histogram.
    // ///
    // /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    // Profile2DPtr bookProfile2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
    //                            const std::string& title="",
    //                            const std::string& xtitle="",
    //                            const std::string& ytitle="",
    //                            const std::string& ztitle="");

    //@}


    /// @name 2D scatter booking
    //@{

    /// @brief Book a 2-dimensional data point set with the given name.
    ///
    /// @note Unlike histogram booking, scatter booking by default makes no
    /// attempt to use reference data to pre-fill the data object. If you want
    /// this, which is sometimes useful e.g. when the x-position is not really
    /// meaningful and can't be extracted from the data, then set the @a
    /// copy_pts parameter to true. This creates points to match the reference
    /// data's x values and errors, but with the y values and errors zeroed...
    /// assuming that there is a reference histo with the same name: if there
    /// isn't, an exception will be thrown.

    Scatter2DPtr & book(Scatter2DPtr & s2d, const string& hname,
                               bool copy_pts=false,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    /// @brief Book a 2-dimensional data point set, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    ///
    /// @note Unlike histogram booking, scatter booking by default makes no
    /// attempt to use reference data to pre-fill the data object. If you want
    /// this, which is sometimes useful e.g. when the x-position is not really
    /// meaningful and can't be extracted from the data, then set the @a
    /// copy_pts parameter to true. This creates points to match the reference
    /// data's x values and errors, but with the y values and errors zeroed.
    Scatter2DPtr & book(Scatter2DPtr & s2d, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                               bool copy_pts=false,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    /// @brief Book a 2-dimensional data point set with equally spaced x-points in a range.
    ///
    /// The y values and errors will be set to 0.
    Scatter2DPtr & book(Scatter2DPtr & s2d, const string& hname,
                               size_t npts, double lower, double upper,
                               const std::string& title="",
                               const std::string& xtitle="",
                               const std::string& ytitle="");

    /// @brief Book a 2-dimensional data point set based on provided contiguous "bin edges".
    ///
    /// The y values and errors will be set to 0.
    Scatter2DPtr & book(Scatter2DPtr & s2d, const string& hname,
                               const std::vector<double>& binedges,
                               const std::string& title,
                               const std::string& xtitle,
                               const std::string& ytitle);

    //@}


  private:
    /// to be used in finalize context only
    class CounterAdapter {
    public:
    	CounterAdapter(double          x) : x_(x                ) {}

	CounterAdapter(const YODA::Counter   & c) : x_(c.val()          ) {}

	// CounterAdapter(CounterPtr     cp) : x_(cp->val()        ) {}

    	CounterAdapter(const YODA::Scatter1D & s) : x_(s.points()[0].x()) {
    		assert( s.numPoints() == 1 || "Can only scale by a single value.");
    	}

    	// CounterAdapter(Scatter1DPtr   sp) : x_(sp->points()[0].x()) {
    	// 	assert( sp->numPoints() == 1 || "Can only scale by a single value.");
    	// }

	operator double() const { return x_; }

    private:
    	double x_;
    };

  public:

    double dbl(double          x) { return x; }
    double dbl(const YODA::Counter   & c) { return c.val(); }
    double dbl(const YODA::Scatter1D & s) {  
    	assert( s.numPoints() == 1 ); 
    	return s.points()[0].x();
    }

    /// @name Analysis object manipulation
    /// @todo Should really be protected: only public to keep BinnedHistogram happy for now...
    //@{

    /// Multiplicatively scale the given counter, @a cnt, by factor @s factor.
    void scale(CounterPtr cnt, CounterAdapter factor);

    /// Multiplicatively scale the given counters, @a cnts, by factor @s factor.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of CounterPtrs
    void scale(const std::vector<CounterPtr>& cnts, CounterAdapter factor) {
      for (auto& c : cnts) scale(c, factor);
    }
    /// @todo YUCK!
    template <std::size_t array_size>
    void scale(const CounterPtr (&cnts)[array_size], CounterAdapter factor) {
      // for (size_t i = 0; i < std::extent<decltype(cnts)>::value; ++i) scale(cnts[i], factor);
      for (auto& c : cnts) scale(c, factor);
    }


    /// Normalize the given histogram, @a histo, to area = @a norm.
    void normalize(Histo1DPtr histo, CounterAdapter norm=1.0, bool includeoverflows=true);

    /// Normalize the given histograms, @a histos, to area = @a norm.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of Histo1DPtrs
    void normalize(const std::vector<Histo1DPtr>& histos, CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& h : histos) normalize(h, norm, includeoverflows);
    }
    /// @todo YUCK!
    template <std::size_t array_size>
    void normalize(const Histo1DPtr (&histos)[array_size], CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& h : histos) normalize(h, norm, includeoverflows);
    }

    /// Multiplicatively scale the given histogram, @a histo, by factor @s factor.
    void scale(Histo1DPtr histo, CounterAdapter factor);

    /// Multiplicatively scale the given histograms, @a histos, by factor @s factor.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of Histo1DPtrs
    void scale(const std::vector<Histo1DPtr>& histos, CounterAdapter factor) {
      for (auto& h : histos) scale(h, factor);
    }
    /// @todo YUCK!
    template <std::size_t array_size>
    void scale(const Histo1DPtr (&histos)[array_size], CounterAdapter factor) {
      for (auto& h : histos) scale(h, factor);
    }


    /// Normalize the given histogram, @a histo, to area = @a norm.
    void normalize(Histo2DPtr histo, CounterAdapter norm=1.0, bool includeoverflows=true);

    /// Normalize the given histograms, @a histos, to area = @a norm.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of Histo2DPtrs
    void normalize(const std::vector<Histo2DPtr>& histos, CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& h : histos) normalize(h, norm, includeoverflows);
    }
    /// @todo YUCK!
    template <std::size_t array_size>
    void normalize(const Histo2DPtr (&histos)[array_size], CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& h : histos) normalize(h, norm, includeoverflows);
    }

    /// Multiplicatively scale the given histogram, @a histo, by factor @s factor.
    void scale(Histo2DPtr histo, CounterAdapter factor);

    /// Multiplicatively scale the given histograms, @a histos, by factor @s factor.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of Histo2DPtrs
    void scale(const std::vector<Histo2DPtr>& histos, CounterAdapter factor) {
      for (auto& h : histos) scale(h, factor);
    }
    /// @todo YUCK!
    template <std::size_t array_size>
    void scale(const Histo2DPtr (&histos)[array_size], CounterAdapter factor) {
      for (auto& h : histos) scale(h, factor);
    }


    /// Helper for counter division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(CounterPtr c1, CounterPtr c2, Scatter1DPtr s) const;

    /// Helper for histogram division with raw YODA objects.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(const YODA::Counter& c1, const YODA::Counter& c2, Scatter1DPtr s) const;


    /// Helper for histogram division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const;

    /// Helper for histogram division with raw YODA objects.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(const YODA::Histo1D& h1, const YODA::Histo1D& h2, Scatter2DPtr s) const;


    /// Helper for profile histogram division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(Profile1DPtr p1, Profile1DPtr p2, Scatter2DPtr s) const;

    /// Helper for profile histogram division with raw YODA objects.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(const YODA::Profile1D& p1, const YODA::Profile1D& p2, Scatter2DPtr s) const;


    /// Helper for 2D histogram division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(Histo2DPtr h1, Histo2DPtr h2, Scatter3DPtr s) const;

    /// Helper for 2D histogram division with raw YODA objects.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(const YODA::Histo2D& h1, const YODA::Histo2D& h2, Scatter3DPtr s) const;


    /// Helper for 2D profile histogram division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(Profile2DPtr p1, Profile2DPtr p2, Scatter3DPtr s) const;

    /// Helper for 2D profile histogram division with raw YODA objects
    ///
    /// @note Assigns to the (already registered) output scatter, @a s.  Preserves the path information of the target.
    void divide(const YODA::Profile2D& p1, const YODA::Profile2D& p2, Scatter3DPtr s) const;


    /// Helper for histogram efficiency calculation.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void efficiency(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const;

    /// Helper for histogram efficiency calculation.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void efficiency(const YODA::Histo1D& h1, const YODA::Histo1D& h2, Scatter2DPtr s) const;


    /// Helper for histogram asymmetry calculation.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void asymm(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const;

    /// Helper for histogram asymmetry calculation.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void asymm(const YODA::Histo1D& h1, const YODA::Histo1D& h2, Scatter2DPtr s) const;


    /// Helper for converting a differential histo to an integral one.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void integrate(Histo1DPtr h, Scatter2DPtr s) const;

    /// Helper for converting a differential histo to an integral one.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void integrate(const Histo1D& h, Scatter2DPtr s) const;

    //@}


  public:

    /// List of registered analysis data objects
    const vector<MultiweightAOPtr>& analysisObjects() const {
      return _analysisobjects;
    }


  protected:

    /// @name Data object registration, and removal
    //@{

    /// Register a data object in the histogram system
    void addAnalysisObject(const MultiweightAOPtr & ao);

    /// Unregister a data object from the histogram system (by name)
    void removeAnalysisObject(const std::string& path);

    /// Unregister a data object from the histogram system (by pointer)
    void removeAnalysisObject(const MultiweightAOPtr & ao);

    //@}


  private:

    /// Name passed to constructor (used to find .info analysis data file, and as a fallback)
    string _defaultname;

    /// Pointer to analysis metadata object
    unique_ptr<AnalysisInfo> _info;

    /// Storage of all plot objects
    /// @todo Make this a map for fast lookup by path?
    vector<MultiweightAOPtr> _analysisobjects;

    /// @name Cross-section variables
    //@{
    double _crossSection;
    //@}

    /// The controlling AnalysisHandler object.
    AnalysisHandler* _analysishandler;

    /// Collection of cached refdata to speed up many autobookings: the
    /// reference data file should only be read once.
    mutable std::map<std::string, YODA::AnalysisObjectPtr> _refdata;


  private:

    /// @name Utility functions
    //@{

    /// Get the reference data for this paper and cache it.
    void _cacheRefData() const;

    //@}


    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    Analysis& operator=(const Analysis&);

  };


}


// Include definition of analysis plugin system so that analyses automatically see it when including Analysis.hh
#include "Rivet/AnalysisBuilder.hh"

/// @def DECLARE_RIVET_PLUGIN
/// Preprocessor define to prettify the global-object plugin hook mechanism.
#define DECLARE_RIVET_PLUGIN(clsname) Rivet::AnalysisBuilder<clsname> plugin_ ## clsname

/// @def DECLARE_ALIASED_RIVET_PLUGIN
/// Preprocessor define to prettify the global-object plugin hook mechanism, with an extra alias name for this analysis.
// #define DECLARE_ALIASED_RIVET_PLUGIN(clsname, alias) Rivet::AnalysisBuilder<clsname> plugin_ ## clsname ## ( ## #alias ## )
#define DECLARE_ALIASED_RIVET_PLUGIN(clsname, alias) DECLARE_RIVET_PLUGIN(clsname)( #alias )

/// @def DEFAULT_RIVET_ANALYSIS_CONSTRUCTOR
/// Preprocessor define to prettify the manky constructor with name string argument
#define DEFAULT_RIVET_ANALYSIS_CONSTRUCTOR(clsname) clsname() : Analysis(# clsname) {}

/// @def DEFAULT_RIVET_ANALYSIS_CTOR
/// Slight abbreviation for DEFAULT_RIVET_ANALYSIS_CONSTRUCTOR
#define DEFAULT_RIVET_ANALYSIS_CTOR(clsname) DEFAULT_RIVET_ANALYSIS_CONSTRUCTOR(clsname)



#endif
