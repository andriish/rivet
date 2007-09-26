// -*- C++ -*-
#ifndef RIVET_AnalysisLoader_HH
#define RIVET_AnalysisLoader_HH

#include "Rivet/Rivet.hh"

// Forward declaration
namespace Rivet {  
  class Analysis;
}

// Typedefs for dlopen() magic
#include <map>
#include <string>
typedef Rivet::Analysis* (*anabuilder_fn)();
typedef std::map<std::string,anabuilder_fn> AnalysisBuilders;
typedef AnalysisBuilders (*anabuilders_fn)();


namespace Rivet {  

  class AnalysisLoader {
  public:
    
    /// Get all the available analyses' names.
    static set<string> getAllAnalysisNames();

    /// Get an analysis by name. The returned Analysis object will be obtained 
    /// from the Analysis::create() factory method via some dlopen() magic.
    /// Warning: a name arg which matches no known analysis will return a null
    /// pointer. Check your return values before using them!
    static Analysis* getAnalysis(const string& analysisname);

    /// Get all the available analyses.
    static set<Analysis*> getAllAnalyses();

    /// Load the analysis builder functions from a named shared object file.
    static AnalysisBuilders& loadAnalysisBuildersFromFile(const string& filename, AnalysisBuilders& builders);

    /// Close the dlopen()ed libraries.    
    static void closeAnalysisBuilders();

    /// Load the analysis builder functions from shared object files in a named directory.
    static AnalysisBuilders& loadAnalysisBuildersFromDir(const string& dirname, AnalysisBuilders& builders); 

    /// Load the analysis builder functions from shared object files in named directories. 
    static AnalysisBuilders& loadAnalysisBuildersFromDirs(const vector<string>& dirnames, AnalysisBuilders& builders);

  private:
    /// Load the available analyses at runtime.
    static void loadAnalyses();
    
    /// Static flag for determining when to run the load method.
    static bool _loaded;
    
    /// The named factory functions of all available analyses.
    static AnalysisBuilders _analysisbuilders;

    /// The named factory functions of all available analyses.
    static set<void*> _handles;
    
  };
  
}

#endif
