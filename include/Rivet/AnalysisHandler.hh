// -*- C++ -*-
#ifndef RIVET_RivetHandler_HH
#define RIVET_RivetHandler_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/AnalysisHandler.fhh"
#include "Rivet/Event.hh"
#include "Rivet/Analysis.hh"
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
    /// @param basefilename the name of the file (no extension) where histograms are to be stored.
    /// @param storetype a string indicating to the AIDA analysis factory
    ///   how to store the histograms. Which strings are allowed depends on
    ///   actual AIDA implementation used. To output in standard AIDA XML the
    ///   string is typically "xml".
    /// @param afac an AIDA analysis factory object. The caller must make
    ///   sure that the lifetime of the factory object exceeds the AnalysisHandler
    ///   object.
    AnalysisHandler(AIDA::IAnalysisFactory& afac, string basefilename="Rivet", HistoFormat storetype=AIDAML);

    /// Make a Rivet handler with a set base filename and store type.
    AnalysisHandler(string basefilename="Rivet", HistoFormat storetype=AIDAML);

    /// The destructor is not virtual as this class should not be inherited from.
    ~AnalysisHandler() { }
    //@}

  private:

    /// Do the initialisation of the AIDA analysis factories.
    void setupFactories(string basefilename, HistoFormat storetype);

    /// Convert any IHistogram1D objects in the AIDA tree to IDataPointSet objects.
    void normalizeTree(AIDA::ITree& tree);

    /// Get a logger object.
    Log& getLog();

  public:

    /// Get the number of events seen. Should only really be used by external
    /// steering code or analyses in the finalize phase.
    size_t numEvents() const { return _numEvents; }

    /// Get the sum of the event weights seen - the weighted equivalent of the
    /// number of events. Should only really be used by external steering code
    /// or analyses in the finalize phase.
    double sumOfWeights() const { return _sumOfWeights; }

    /// Add an analysis to the run list using its name. The actual Analysis 
    /// to be used will be obtained via AnalysisHandler::getAnalysis(string).
    /// If no matching analysis is found, no analysis is added (i.e. the
    /// null pointer is checked and discarded.
    AnalysisHandler& addAnalysis(const string& analysisname);

    /// Add analyses to the run list using their names. The actual {@link
    /// Analysis}' to be used will be obtained via
    /// AnalysisHandler::addAnalysis(string), which in turn uses
    /// AnalysisHandler::getAnalysis(string). If no matching analysis is found
    /// for a given name, no analysis is added, but also no error is thrown.
    AnalysisHandler& addAnalyses(const vector<string>& analysisnames) {
      for (vector<string>::const_iterator aname = analysisnames.begin();
           aname != analysisnames.end(); ++aname) {
        addAnalysis(*aname);
      }
      return *this;
    }

    /// Add an analysis to the run list by supplying a "template" analysis.
    /// @todo Is there a good reason to not allow "direct" submission?
    template <typename A>
    AnalysisHandler& addAnalysis(const A& analysis) {
      A* a = new A(analysis);
      a->_theHandler = this;
      _analyses.insert(a);
      return *this;
    }

    /// Initialize a run. If this run is to be joined together with other
    /// runs, \a N should be set to the total number of runs to be
    /// combined, and \a i should be the index of this run. This function
    /// will initialize the histogram factory and then call the
    /// AnalysisBase::init() function of all included analysis objects.
    void init(int i=0, int N=0);


    /// Analyze the given \a event. This function will call the
    /// AnalysisBase::analyze() function of all included analysis objects.
    void analyze(const GenEvent& event);


    /// Finalize a run. This function first calls the AnalysisBase::finalize()
    /// functions of all included analysis objects and converts all histograms
    /// to AIDA DataPointSet objects in the AIDA tree. Using the histogram tree
    /// for further analysis or writing to file is left to the API user.
    void finalize();


    /// The AIDA analysis factory.
    AIDA::IAnalysisFactory& analysisFactory() {
      return *_theAnalysisFactory;
    }
    

    /// The AIDA tree object.
    AIDA::ITree& tree() {
      return *_theTree;
    }

    
    /// The AIDA histogram factory.
    AIDA::IHistogramFactory& histogramFactory() {
      return *_theHistogramFactory;
    }


    /// The AIDA histogram factory.
    AIDA::IDataPointSetFactory& datapointsetFactory() {
      return *_theDataPointSetFactory;
    }


    /// Is cross-section information required by at least one child analysis?
    bool needCrossSection() {
      bool rtn = false;
      foreach (Analysis* a, _analyses) {
        if (!rtn) rtn = a->needsCrossSection();
        if (rtn) break;
      }
      return rtn;
    }


    /// Set the cross-section for the process being generated.    
    AnalysisHandler& setCrossSection(const double & xs) {
      foreach (Analysis* a, _analyses) {
        a->setCrossSection(xs);
      }
      return *this;
    }


  private:

    /// The collection of Analysis objects to be used.
    set<Analysis*> _analyses;
    
    /// If non-zero the number of runs to be combined into one analysis.
    int _nRun;

    /// If non-zero, the index of this run, if a part of several runs to
    /// be combined into one.
    int _iRun;

    /// Number of events seen.
    size_t _numEvents;

    /// Sum of event weights seen.
    double _sumOfWeights;

    /// The AIDA analysis factory.
    AIDA::IAnalysisFactory* _theAnalysisFactory;

    /// The AIDA tree object.
    AIDA::ITree* _theTree;

    /// The AIDA histogram factory.
    AIDA::IHistogramFactory* _theHistogramFactory;

    /// The AIDA data point set factory.
    AIDA::IDataPointSetFactory* _theDataPointSetFactory;

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
