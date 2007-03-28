// -*- C++ -*-
#ifndef RIVET_RivetHandler_H
#define RIVET_RivetHandler_H
//
// This is the declaration of the RivetHandler class.
//

#include "Rivet/Rivet.hh"
#include "Rivet/RivetHandler.fhh"
#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Tools/Event/Event.hh"


namespace AIDA {
  class IAnalysisFactory;
  class IHistogramFactory;
  class ITree;
}


namespace Rivet {

  /**
   * RivetHandler is a concrete class which handles a number of analysis
   * objects to be applied to generated events.
   */
  class RivetHandler {

  public:

    /** Typedef a vector of analysis objects. */
    typedef vector<AnalysisPtr> AnalysisVector;

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * The standard constructor.
     * @param filename the name of the file (no extension) where histograms are to be stored.
     * @param storetype a string indicating to the AIDA analysis factory
     * how to store the histograms. Which strings are allowed depends on
     * actual AIDA implementation used. To output in standard AIDA XML the
     * string is typically "xml".
     * @param afac an AIDA analysis factory object. The caller must make
     * sure that the lifetime of the factory object exceeds the RivetHandler
     * object.
     */
    RivetHandler(AIDA::IAnalysisFactory& afac, std::string basefilename="Rivet", HistoFormat storetype=AIDAML);

    /// Make a Rivet handler with a set base filename and store type.
    /// @todo storetype -> enum
    RivetHandler(std::string basefilename="Rivet", HistoFormat storetype=AIDAML);

    /**
     * The destructor is not virtual as this class should not be
     * inherited from.
     */
    ~RivetHandler();
    //@}

  private:

    void setupFactories(string basefilename, HistoFormat storetype);

  public:


    /// Add an object of base class Analysis to the list of analysis
    /// objects to be used in a run. Note that the object will be copied,
    /// and the argument object given will not be used and can be
    /// discarded directly.
    template <typename A>
    inline RivetHandler& addAnalysis(const A& analysis) {
      analysisVector.push_back(new A(analysis));
      analysisVector.back()->theHandler = this;
      return *this;
    }


    /// Add a collection of analyses
    /// @todo Make the addAnalyses(vector<Analysis*>) method work
    inline RivetHandler& addAnalyses(const std::vector<Analysis*>& analyses) { 
      throw runtime_error("addAnalyses(vector<Analysis*>) doesn't yet work...");
      //for (std::vector<Analysis*>::const_iterator ana = analyses.begin();
      //     ana != analyses.end(); ++ana) {
      //  Analysis& atemp = **ana;
      //  addAnalysis(atemp);
      //}
      return *this;
    }


    /// Add an object of base class Analysis to the list of analysis objects to
    /// be used in a run. The Analysis will be obtained from the
    /// Analysis::getAnalysis() factory method, according to the argument enum.
    inline RivetHandler& addAnalysis(const AnalysisName analysisname) { 
      Analysis& analysis = Analysis::getAnalysis(analysisname);
      analysisVector.push_back(&analysis);
      analysisVector.back()->theHandler = this;
      return *this;
    }


    /// Add a collection of analyses via their analysis name enums
    inline RivetHandler& addAnalyses(cAnalysisList& analysisnames) { 
      for (cAnalysisList::const_iterator ananame = analysisnames.begin();
           ananame != analysisnames.end(); ++ananame) {
        addAnalysis(*ananame);
      }
      return *this;
    }

    
    /// Initialize a run. If this run is to be joined together with other
    /// runs, \a N should be set to the total number of runs to be
    /// combined, and \a i should be the index of this run. This function
    /// will initialize the histogram factory and then call the
    /// AnalysisBase::init() function of all included analysis objects.
    void init(int i = 0, int N = 0);


    /// Analyze the given \a event. This function will call the
    /// AnalysisBase::analyze() function of all included analysis objects.
    void analyze(const GenEvent & event);


    /// Finalize a run. This function first calls the
    /// AnalysisBase::finalize() functions of all included analysis
    /// objects and then writes out all histograms to a file.
    void finalize();


    /// Return a RivetInfo object containing all parameters of all
    /// included analysis handlers and all their included Projections.
    RivetInfo info() const;


    /// The AIDA analysis factory.
    inline AIDA::IAnalysisFactory& analysisFactory() {
      return *theAnalysisFactory;
    }
    

    /// The AIDA tree object.
    inline AIDA::ITree& tree() {
      return *theTree;
    }
    
    
    /// The AIDA histogram factory.
    inline AIDA::IHistogramFactory& histogramFactory() {
      return *theHistogramFactory;
    }
    


  private:

    /// The vector of Analysis objects to be used.
    AnalysisVector analysisVector;

    /// If non-zero the number of runs to be combined into one analysis.
    int nRun;

    /// If non-zero, the index of this run, if a part of several runs to
    /// be combined into one.
    int iRun;

    /// The AIDA analysis factory.
    AIDA::IAnalysisFactory* theAnalysisFactory;

    /// The AIDA tree object.
    AIDA::ITree* theTree;

    /// The AIDA histogram factory.
    AIDA::IHistogramFactory* theHistogramFactory;

  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    RivetHandler & operator=(const RivetHandler &);

    /// The copy constructor is private and must never be called.  In
    /// fact, it should not even be implemented.
    inline RivetHandler(const RivetHandler &);

  };


}

#endif /* RIVET_RivetHandler_H */
