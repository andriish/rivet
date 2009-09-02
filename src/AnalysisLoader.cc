// -*- C++ -*-
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/osdir.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis.hh"
#include <dlfcn.h>

namespace Rivet {


  // Initialise static ptr collection
  AnalysisLoader::AnalysisBuilderMap AnalysisLoader::_ptrs;


  vector<string> AnalysisLoader::analysisNames() {
    _loadAnalysisPlugins();
    vector<string> names;
    foreach (const AnalysisBuilderMap::value_type& p, _ptrs) names += p.first;
    return names;
  }


  Analysis* AnalysisLoader::getAnalysis(const string& analysisname) { 
    _loadAnalysisPlugins();
    AnalysisBuilderMap::const_iterator ai = _ptrs.find(analysisname);
    if (ai == _ptrs.end()) return 0;
    return ai->second->mkAnalysis();
  }


  vector<Analysis*> AnalysisLoader::getAllAnalyses() {
    _loadAnalysisPlugins();
    vector<Analysis*> analyses;
    foreach (const AnalysisBuilderMap::value_type& p, _ptrs) {
      analyses += p.second->mkAnalysis();
    }
    return analyses;
  }


  void AnalysisLoader::_registerBuilder(const AnalysisBuilderBase* ab) {
    if (!ab) return;
    const string name = ab->name();
    if (_ptrs.find(name) != _ptrs.end()) {
      // Protect against overwriting analyses by throwing an error
      /// @todo Tidy this up!
      Log::getLog("Rivet.AnalysisLoader") << Log::ERROR << "Tried to register a second plugin analysis called '" << name << "'" << endl;
      exit(1);
    }
    Log::getLog("Rivet.AnalysisLoader") << Log::TRACE << "Registering a plugin analysis called '" << name << "'" << endl;
    _ptrs[name] = ab;
  }
  

  void AnalysisLoader::_loadAnalysisPlugins() {
    // Only run once
    if (!_ptrs.empty()) return;

    #ifdef LIB_SUFFIX
    const string libsuffix(LIB_SUFFIX);
    #else
    const string libsuffix(".so");
    #endif

    // Build the list of directories to search
    vector<string> dirs;
    char* env = 0;
    // First use the Rivet analysis path variable
    env = getenv("RIVET_ANALYSIS_PATH");
    if (env) dirs += split(env);
    // Then the Rivet library install path
    dirs += getLibPath();    
    // And then the user's (non-system) library path
    env = getenv("LD_LIBRARY_PATH");
    if (env) dirs += split(env);
    // Use the current dir, and its libtool subdir, too.
    dirs += ".";
    dirs += ".libs";

    // Find plugin module library files
    vector<string> pluginfiles;
    foreach (const string& d, dirs) {
      if (d.empty()) continue;
      oslink::directory dir(d);
      while (dir) {
        string filename = dir.next();
        // Require that name *starts* with 'Rivet' with new loader
        if (filename.find("Rivet") != 0) continue;
        size_t posn = filename.find(libsuffix);
        if (posn == string::npos || posn != filename.length()-libsuffix.length()) continue;
        /// @todo Make sure this is an abs path
        /// @todo Sys-dependent path separator instead of "/"
        const string path = d + "/" + filename;
        // Ensure no duplicate paths
        if (find(pluginfiles.begin(), pluginfiles.end(), path) == pluginfiles.end()) {
          pluginfiles += path;
        }
      }
    }

    // Load the plugin files
    Log::getLog("Rivet.AnalysisLoader") << Log::TRACE << "Candidate analysis plugin libs: " << pluginfiles << endl;
    foreach (const string& pf, pluginfiles) {
      Log::getLog("Rivet.AnalysisLoader") << Log::TRACE << "Trying to load plugin analyses from file " << pf << endl;
      /// @todo Still needs lazy loading with new plugin design?
      void* handle = dlopen(pf.c_str(), RTLD_LAZY);
      if (!handle) {
        Log::getLog("Rivet.AnalysisLoader") << Log::WARN << "Cannot open " << pf << ": " << dlerror() << endl;
        continue;
      }
    }
  }


}
