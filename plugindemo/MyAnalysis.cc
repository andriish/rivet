#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {
  
  class MyAnalysis : public Analysis {
    
  public:
    
    /// Default constructor
    MyAnalysis() {
      const FinalState cnfs;
      addProjection(cnfs, "CNFS");
      const ChargedFinalState cfs;
      addProjection(cfs, "CFS");
    }
    
    
    /// Factory method
    static Analysis* create() { 
      return new MyAnalysis(); 
    }
    
    
    /// Return the name of this analysis
    string name() const {
      return "MyAnalysis";
    }
    
    /// Get the SPIRES ID code
    string spiresId() const {
      return "NONE";
    }
    
    /// Get a description of the analysis.
    string description() const {
      return "A do-nothing analysis for demonstrating how to make a plugin";
    }
    
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "NONE";
    }
    
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "NONE";
    }
    
    
    /// @name Analysis methods
    //@{
    void init() { 
      // No histos, so nothing to do!
    }
    
    void analyze(const Event& event) {
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      const FinalState& cnfs = applyProjection<FinalState>(event, "CNFS");
      getLog() << Log::DEBUG << "Total multiplicity            = " 
               << cnfs.particles().size() << endl;
      getLog() << Log::DEBUG << "Total charged multiplicity    = " 
               << cfs.particles().size()  << endl;
    }
    
    void finalize() {
      // No histos, so nothing to do!
    }
    //@}
    
  };
  
  
  extern "C" {
    AnalysisBuilders getAnalysisBuilders() {
      AnalysisBuilders fns;
      fns["MYANALYSIS"] = Rivet::MyAnalysis::create;
      return fns;
    }
  }

  
}
