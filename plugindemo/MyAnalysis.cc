#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {
  

  class MyAnalysis : public Analysis {
    
  public:
    
    /// Default constructor
    MyAnalysis() : Analysis("MYANALYSIS") {
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
    
    /// Summary of analysis
    string summary() const{
      return "An example plugin analysis. Prints the total and charged multiplicities to the screen for each event.";
    }
    
    /// Beam conditions for this analysis
    string runInfo() const{
      return "Anything!";
    }
    
    string collider() const{
      return "Any";
    }
    
    vector<string> authors() const{
      return vector<string>();
    }
    
    vector<string> references() const{
      return vector<string>();
    }
    
    /// @name Analysis methods
    //@{
    void init() { 
      const FinalState cnfs;
      addProjection(cnfs, "CNFS");
      const ChargedFinalState cfs;
      addProjection(cfs, "CFS");
    }
    
    void analyze(const Event& event) {
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      const FinalState& cnfs = applyProjection<FinalState>(event, "CNFS");
      getLog() << Log::DEBUG << "Total multiplicity = " << cnfs.size() << endl;
      getLog() << Log::DEBUG << "Total charged multiplicity = " << cfs.size()  << endl;
    }
    
    void finalize() {
      // No histos, so nothing to do!
    }
    //@}
    
  };
  
  

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MyAnalysis> plugin_MyAnalysis;
  
}
