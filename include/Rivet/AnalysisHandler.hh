// -*- C++ -*-

#ifndef RIVET_RivetHandler_H
#define RIVET_RivetHandler_H

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/RivetAIDA.fhh"
#include "Rivet/AnalysisHandler.fhh"
#include "Rivet/Event.hh"
#include "Rivet/Analysis/Analysis.hh"


namespace Rivet {


  /// A class which handles a number of analysis objects to be applied to 
  /// generated events. An {@link Analysis}' AnalysisHandler is also responsible
  /// for handling the final writing-out of histograms.
  class AnalysisHandler {

  public:

    /// Typedef a vector of analysis objects.
    typedef vector<AnalysisPtr> AnalysisVector;

  public:

    /// @name Standard constructors and destructors. */
    //@{
    /// The standard constructor.
    /// @param filename the name of the file (no extension) where histograms are to be stored.
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
    ~AnalysisHandler();
    //@}

  private:

    /// Do the initialisation of the AIDA analysis factories.
    void setupFactories(string basefilename, HistoFormat storetype);

    /// Convert any IHistogram1D objects in the AIDA tree to IDataPointSet objects.
    void normalizeTree(AIDA::ITree& tree);

    /// Get a logger object.
    Log& getLog();

  public:

    /// Add an object of base class Analysis to the list of analysis
    /// objects to be used in a run. Note that the object will be copied,
    /// and the argument object given will not be used and can be
    /// discarded directly.
    template <typename A>
    inline AnalysisHandler& addAnalysis(const A& analysis) {
      /// @todo Check that there aren't any duplicate analyses?
      _analysisVector.push_back(new A(analysis));
      _analysisVector.back()->_theHandler = this;
      return *this;
    }


    /// Add a collection of analyses. Note that this method doesn't make a 
    /// copy of the passed analyses, since it can't dereference the abstract 
    /// Analysis* pointer - be careful not to use the parameter analyses 
    /// afterwards - Rivet will have taken ownership of them.
    inline AnalysisHandler& addAnalyses(const vector<Analysis*>& analyses) { 
      throw runtime_error("addAnalyses(vector<Analysis*>) doesn't yet work...");
      for (vector<Analysis*>::const_iterator ana = analyses.begin(); ana != analyses.end(); ++ana) {
        _analysisVector.push_back(*ana);
        _analysisVector.back()->_theHandler = this;
      }
      return *this;
    }


    /// Add an object of base class Analysis to the list of analysis objects to
    /// be used in a run. The Analysis will be obtained from the
    /// Analysis::getAnalysis() factory method, according to the argument enum.
    inline AnalysisHandler& addAnalysis(const AnalysisName analysisname) { 
      Analysis& analysis = Analysis::getAnalysis(analysisname);
      _analysisVector.push_back(&analysis);
      _analysisVector.back()->_theHandler = this;
      return *this;
    }


    /// Add a collection of analyses via their analysis name enums
    inline AnalysisHandler& addAnalyses(cAnalysisList& analysisnames) { 
      for (cAnalysisList::const_iterator ananame = 
             analysisnames.begin(); ananame != analysisnames.end(); ++ananame) {
        addAnalysis(*ananame);
      }
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


    /// Finalize a run. This function first calls the
    /// AnalysisBase::finalize() functions of all included analysis
    /// objects and then writes out all histograms to a file.
    void finalize();


    /// The AIDA analysis factory.
    inline AIDA::IAnalysisFactory& analysisFactory() {
      return *_theAnalysisFactory;
    }
    

    /// The AIDA tree object.
    inline AIDA::ITree& tree() {
      return *_theTree;
    }
    
    
    /// The AIDA histogram factory.
    inline AIDA::IHistogramFactory& histogramFactory() {
      return *_theHistogramFactory;
    }
    

    /// The AIDA histogram factory.
    inline AIDA::IDataPointSetFactory& datapointsetFactory() {
      return *_theDataPointSetFactory;
    }



  private:

    /// The vector of Analysis objects to be used.
    AnalysisVector _analysisVector;

    /// If non-zero the number of runs to be combined into one analysis.
    int _nRun;

    /// If non-zero, the index of this run, if a part of several runs to
    /// be combined into one.
    int _iRun;

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
