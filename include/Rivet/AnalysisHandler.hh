// -*- C++ -*-
#ifndef RIVET_RivetHandler_HH
#define RIVET_RivetHandler_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet {


  // Forward declaration and smart pointer for Analysis
  class Analysis;
  typedef std::shared_ptr<Analysis> AnaHandle;


  // Needed to make smart pointers compare equivalent in the STL set
  struct CmpAnaHandle {
    bool operator() (const AnaHandle& a, const AnaHandle& b) const {
      return a.get() < b.get();
    }
  };


  /// A class which handles a number of analysis objects to be applied to
  /// generated events. An {@link Analysis}' AnalysisHandler is also responsible
  /// for handling the final writing-out of histograms.
  class AnalysisHandler {
  public:

    /// @name Constructors and destructors. */
    //@{

    /// Preferred constructor, with optional run name.
    AnalysisHandler(const string& runname="");

    /// @brief Destructor
    /// The destructor is not virtual, as this class should not be inherited from.
    ~AnalysisHandler();

    //@}


  private:

    /// Get a logger object.
    Log& getLog() const;


  public:

    /// @name Run properties
    //@{

    /// Get the name of this run.
    string runName() const;


    /// Get the number of events seen. Should only really be used by external
    /// steering code or analyses in the finalize phase.
    size_t numEvents() const;

    /// Get the sum of the event weights seen - the weighted equivalent of the
    /// number of events. Should only really be used by external steering code
    /// or analyses in the finalize phase.
    double sumOfWeights() const {
        return _eventCounter->sumW();
    }

    const vector<string>& weightNames() const {
        return _weightNames;
    }

    size_t numWeights() const {
        return _weightNames.size();
    }

    /// Set the cross-section for the process being generated.
    AnalysisHandler& setCrossSection(double xs, double xserr);

    /// Get the cross-section known to the handler.
    Scatter1DPtr crossSection() const {
      return _xs;
    }

    double nominalCrossSection() const {
      _xs.get()->setActiveWeightIdx(_defaultWeightIdx);
      const YODA::Scatter1D::Points& ps = _xs->points();
      if (ps.size() != 1) {
        string errMsg = "cross section missing when requesting nominal cross section";
        throw Error(errMsg);
      }
      double xs = ps[0].x();
      _xs.get()->unsetActiveWeight();
      return xs;
    }


    /// Set the beam particles for this run
    AnalysisHandler& setRunBeams(const ParticlePair& beams) {
      _beams = beams;
      MSG_DEBUG("Setting run beams = " << beams << " @ " << sqrtS()/GeV << " GeV");
      return *this;
    }

    /// Get the beam particles for this run, usually determined from the first event.
    const ParticlePair& beams() const { return _beams; }

    /// Get beam IDs for this run, usually determined from the first event.
    /// @deprecated Use standalone beamIds(ah.beams()), to clean AH interface
    PdgIdPair beamIds() const;

    /// Get energy for this run, usually determined from the first event.
    /// @deprecated Use standalone sqrtS(ah.beams()), to clean AH interface
    double sqrtS() const;

    /// Setter for _ignoreBeams
    void setIgnoreBeams(bool ignore=true);

    /// Setter for _skipWeights
    void skipMultiWeights(bool ignore=false);

    //@}


    /// @name Handle analyses
    //@{

    /// Get a list of the currently registered analyses' names.
    std::vector<std::string> analysisNames() const;

    /// Get the collection of currently registered analyses.
    const std::set<AnaHandle, CmpAnaHandle>& analyses() const {
      return _analyses;
    }

    /// Get a registered analysis by name.
    const AnaHandle analysis(const std::string& analysisname) const;


    /// Add an analysis to the run list by object
    AnalysisHandler& addAnalysis(Analysis* analysis);

    /// @brief Add an analysis to the run list using its name.
    ///
    /// The actual Analysis to be used will be obtained via
    /// AnalysisLoader::getAnalysis(string).  If no matching analysis is found,
    /// no analysis is added (i.e. the null pointer is checked and discarded.
    AnalysisHandler& addAnalysis(const std::string& analysisname);

    /// @brief Add analyses to the run list using their names.
    ///
    /// The actual {@link Analysis}' to be used will be obtained via
    /// AnalysisHandler::addAnalysis(string), which in turn uses
    /// AnalysisLoader::getAnalysis(string). If no matching analysis is found
    /// for a given name, no analysis is added, but also no error is thrown.
    AnalysisHandler& addAnalyses(const std::vector<std::string>& analysisnames);


    /// Remove an analysis from the run list using its name.
    AnalysisHandler& removeAnalysis(const std::string& analysisname);

    /// Remove analyses from the run list using their names.
    AnalysisHandler& removeAnalyses(const std::vector<std::string>& analysisnames);

    //@}


    /// @name Main init/execute/finalise
    //@{

    /// Initialize a run, with the run beams taken from the example event.
    void init(const GenEvent& event);

    /// @brief Analyze the given \a event by reference.
    ///
    /// This function will call the AnalysisBase::analyze() function of all
    /// included analysis objects.
    void analyze(const GenEvent& event);

    /// @brief Analyze the given \a event by pointer.
    ///
    /// This function will call the AnalysisBase::analyze() function of all
    /// included analysis objects, after checking the event pointer validity.
    void analyze(const GenEvent* event);

    /// Finalize a run. This function calls the AnalysisBase::finalize()
    /// functions of all included analysis objects.
    void finalize();

    //@}


    /// @name Histogram / data object access
    //@{

    /// Add a vector of analysis objects to the current state.
    void addData(const std::vector<YODA::AnalysisObjectPtr>& aos);

    /// Read analysis plots into the histo collection (via addData) from the named file.
    void readData(const std::string& filename);

    /// Get all analyses' plots as a vector of analysis objects.
    std::vector<YODA::AnalysisObjectPtr> getData() const;
    std::vector<MultiweightAOPtr>        getRivetAOs() const;
    std::vector<YODA::AnalysisObjectPtr> getYodaAOs() const;

    /// Get all analyses' plots as a vector of analysis objects.
    void setWeightNames(const GenEvent& ge);

    /// Do we have named weights?
    bool haveNamedWeights();


    /// Write all analyses' plots (via getData) to the named file.
    void writeData(const std::string& filename) const;

    //@}

    /// Indicate which Rivet stage we're in.
    /// At the moment, only INIT is used to enable booking.
    enum class Stage { OTHER, INIT };

    /// Which stage are we in?
    Stage stage() const { return _stage; }

  private:
    /// Current handler stage
    Stage _stage = Stage::OTHER;

    /// The collection of Analysis objects to be used.
    set<AnaHandle, CmpAnaHandle> _analyses;

    /// @name Run properties
    //@{

    /// Weight names
    std::vector<std::string> _weightNames;
    std::vector<std::valarray<double> > _subEventWeights;
    size_t _numWeightTypes; // always == WeightVector.size()

    /// Run name
    std::string _runname;

    mutable CounterPtr _eventCounter;

    /// Cross-section known to AH
    Scatter1DPtr _xs;

    /// Beams used by this run.
    ParticlePair _beams;

    /// Flag to check if init has been called
    bool _initialised;

    /// Flag whether input event beams should be ignored in compatibility check
    bool _ignoreBeams;

    /// Flag to check if multiweights should be included
    bool _skipWeights;

    /// Current event number
    int _eventNumber;

    size_t _defaultWeightIdx;
    //@}


  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    AnalysisHandler& operator=(const AnalysisHandler&);

    /// The copy constructor is private and must never be called.  In
    /// fact, it should not even be implemented.
    AnalysisHandler(const AnalysisHandler&);

  };


}

#endif
