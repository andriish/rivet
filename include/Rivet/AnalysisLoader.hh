// -*- C++ -*-
#ifndef RIVET_AnalysisLoader_HH
#define RIVET_AnalysisLoader_HH

#include "Rivet/Config/RivetCommon.hh"
#include <map>
#include <string>

namespace Rivet {


  // Forward declarations
  class Analysis;
  class AnalysisBuilderBase;
  class Log;


  /// Internal class which loads and registers analyses from plugin libs
  class AnalysisLoader {
  public:

    /// Get the available analyses' canonical names
    static vector<string> analysisNames();

    /// Get all the available analyses' names, including aliases
    static vector<string> allAnalysisNames();
    /// @deprecated Use allAnalysisNames()
    static vector<string> getAllAnalysisNames() { return allAnalysisNames(); }

    /// Get the standard analyses' names (from a release-specific list file)
    static vector<string> stdAnalysisNames();

    /// Get the map of analysis alias-names to their canonical equivalents
    static map<string,string> analysisNameAliases();


    /// Get an analysis by name.
    /// Warning: a name arg which matches no known analysis will return a null
    /// pointer. Check your return values before using them!
    static unique_ptr<Analysis> getAnalysis(const string& analysisname);

    /// Get all the available analyses.
    static vector<unique_ptr<Analysis>> getAllAnalyses();


  private:

    /// Allow the analysis builders to call the private _registerBuilder function
    friend class AnalysisBuilderBase;

    /// Register a new analysis builder
    static void _registerBuilder(const AnalysisBuilderBase* ab);

    /// Load the available analyses at runtime.
    static void _loadAnalysisPlugins();

    typedef map<string, const AnalysisBuilderBase*> AnalysisBuilderMap;
    static AnalysisBuilderMap _ptrs;
    static AnalysisBuilderMap _aliasptrs;

  };


}

#endif
