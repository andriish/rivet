// -*- C++ -*-
#ifndef RIVET_RivetHandler_HH
#define RIVET_RivetHandler_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/AnalysisHandler.fhh"
#include "Rivet/Analysis.fhh"
#include "Rivet/Event.fhh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {


  /// A class which handles a number of analysis objects to be applied to
  /// generated events. An {@link Analysis}' AnalysisHandler is also responsible
  /// for handling the final writing-out of histograms.
  class AnalysisHandler {

  public:

    /// @name Standard constructors and destructors. */
    //@{
    /// The standard constructor.
    /// @param basefilename the name of the file (no extension) where histograms
    ///   are to be stored.
    /// @param runname optional name of this run, prepended to AIDA data paths.
    /// @param storetype a string indicating to the AIDA analysis factory
    ///   how to store the histograms. Which strings are allowed depends on
    ///   actual AIDA implementation used. To output in standard AIDA XML the
    ///   string is typically "xml".
    /// @param afac an AIDA analysis factory object. The caller must make
    ///   sure that the lifetime of the factory object exceeds the AnalysisHandler
    ///   object.
    AnalysisHandler(AIDA::IAnalysisFactory& afac, string basefilename="Rivet",
                    string runname="", HistoFormat storetype=AIDAML);

    /// Make a Rivet handler with a set base filename and store type.
    AnalysisHandler(string basefilename="Rivet",
                    string runname="", HistoFormat storetype=AIDAML);

    /// The destructor is not virtual as this class should not be inherited from.
    ~AnalysisHandler();

    //@}


  private:

    /// Do the initialisation of the AIDA analysis factories.
    void _setupFactories(string basefilename, HistoFormat storetype);

    /// Convert any IHistogram1D objects in the AIDA tree to IDataPointSet objects.
    void _normalizeTree(AIDA::ITree& tree);

    /// Get a logger object.
    Log& getLog();


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
    double sumOfWeights() const;

    /// Set sum of weights. This is useful if Rivet is steered externally and
    /// the analyses are run for a sub-contribution of the events
    /// (but of course have to be normalised to the total sum of weights)
    void setSumOfWeights(const double& sum);   


    /// Is cross-section information required by at least one child analysis?
    bool needCrossSection() const;

    /// Set the cross-section for the process being generated.
    AnalysisHandler& setCrossSection(double xs);

    /// Get the cross-section known to the handler.
    double crossSection() const {
      return _xs;
    }

    /// Whether the handler knows about a cross-section.
    bool hasCrossSection() const;


    /// Set beams for this run
    AnalysisHandler& setRunBeams(const ParticlePair& beams) { 
      _beams = beams;
      getLog() << Log::DEBUG << "Setting run beams = " << beams
               << " @ " << sqrtS()/GeV << " GeV" << endl;
      return *this;
    }

    /// Get beam IDs for this run, determined from first event
    const ParticlePair& beams() const { 
      return _beams;
    }

    /// Get beam IDs for this run, determined from first event
    BeamPair beamIds() const;

    /// Get energy for this run, determined from first event
    double sqrtS() const;
    
    //@}


    /// @name Handle analyses
    //@{

    /// Get a list of the currently registered analyses' names.
    std::vector<std::string> analysisNames() const;

    /// Get a list of the currently registered analyses' names.
    const std::set<Analysis*>& analyses() const {
      return _analyses;
    }

    /// Add an analysis to the run list using its name. The actual Analysis
    /// to be used will be obtained via AnalysisHandler::getAnalysis(string).
    /// If no matching analysis is found, no analysis is added (i.e. the
    /// null pointer is checked and discarded.
    AnalysisHandler& addAnalysis(const std::string& analysisname);

    /// Remove an analysis from the run list using its name.
    AnalysisHandler& removeAnalysis(const std::string& analysisname);


    /// Add analyses to the run list using their names. The actual {@link
    /// Analysis}' to be used will be obtained via
    /// AnalysisHandler::addAnalysis(string), which in turn uses
    /// AnalysisHandler::getAnalysis(string). If no matching analysis is found
    /// for a given name, no analysis is added, but also no error is thrown.
    AnalysisHandler& addAnalyses(const std::vector<std::string>& analysisnames);

    /// Remove analyses from the run list using their names.
    AnalysisHandler& removeAnalyses(const std::vector<std::string>& analysisnames);


    /// Add an analysis to the run list by object
    AnalysisHandler& addAnalysis(Analysis* analysis);

    /// Remove beam-incompatible analyses from the run list.
    /// @todo Do this automatically in the init phase (including energies) and deprecate explicit use
    AnalysisHandler& removeIncompatibleAnalyses(const BeamPair& beams);

    //@}


    /// @name Main init/execute/finalise
    //@{

    /// @deprecated Obsolete method, kept only for backwards compatibility
    void init() {}


    /// Initialize a run, with the run beams taken from the example event.
    void init(const GenEvent& event);


    /// Analyze the given \a event. This function will call the
    /// AnalysisBase::analyze() function of all included analysis objects.
    void analyze(const GenEvent& event);


    /// Finalize a run. This function first calls the AnalysisBase::finalize()
    /// functions of all included analysis objects and converts all histograms
    /// to AIDA DataPointSet objects in the AIDA tree. Using the histogram tree
    /// for further analysis or writing to file is left to the API user.
    void finalize();

    //@}


    /// @name AIDA factories etc.
    /// @deprecated All this will be removed when histogramming is overhauled
    //@{

    /// The AIDA analysis factory.
    AIDA::IAnalysisFactory& analysisFactory();


    /// Commit the AIDA tree to file.
    void commitData();
 

    /// The AIDA tree object.
    AIDA::ITree& tree();

 
    /// The AIDA histogram factory.
    AIDA::IHistogramFactory& histogramFactory();


    /// The AIDA histogram factory.
    AIDA::IDataPointSetFactory& datapointsetFactory();

    //@}


  private:

    /// The collection of Analysis objects to be used.
    set<Analysis*> _analyses;


    /// @name Run properties
    //@{

    /// Run name
    std::string _runname;
 
    /// Number of events seen.
    size_t _numEvents;

    /// Sum of event weights seen.
    double _sumOfWeights;

    /// Cross-section known to AH
    double _xs;

    /// Beams used by this run.
    ParticlePair _beams;

    /// Flag to check if init has been called
    bool _initialised;

    //@}


    /// @name AIDA factory handles
    //@{

    /// The AIDA analysis factory.
    AIDA::IAnalysisFactory* _theAnalysisFactory;

    /// The AIDA tree object.
    AIDA::ITree* _theTree;

    /// The AIDA histogram factory.
    AIDA::IHistogramFactory* _theHistogramFactory;

    /// The AIDA data point set factory.
    AIDA::IDataPointSetFactory* _theDataPointSetFactory;

    //@}


  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    AnalysisHandler& operator=(const AnalysisHandler&);

    /// The copy constructor is private and must never be called.  In
    /// fact, it should not even be implemented.
    AnalysisHandler(const AnalysisHandler&);

    static void initializeParticleNames();

  };


}

#endif
