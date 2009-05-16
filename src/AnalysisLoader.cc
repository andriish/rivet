// -*- C++ -*-
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/osdir.hh"
#include <dlfcn.h>

namespace Rivet {


  bool AnalysisLoader::_loaded = false;
  AnalysisBuilders AnalysisLoader::_analysisbuilders;
  set<void*> AnalysisLoader::_handles;


  set<string> AnalysisLoader::getAllAnalysisNames() {
    if (!_loaded) loadAnalyses();
    set<string> names;
    foreach (const AnalysisBuilders::value_type& ab, _analysisbuilders) {
      names.insert(ab.first);
    }
    return names;

  }


  Analysis* AnalysisLoader::getAnalysis(const string& analysisname) { 
    if (!_loaded) loadAnalyses();
    AnalysisBuilders::iterator ab = _analysisbuilders.find(toUpper(analysisname));
    if (ab != _analysisbuilders.end()) {
      return ab->second();
    } else {
      return 0;
    }
  }


  set<Analysis*> AnalysisLoader::getAllAnalyses() {
    if (!_loaded) loadAnalyses();
    set<Analysis*> analyses;
    const set<string> names = getAllAnalysisNames();
    foreach (const string& n, names) {
      analyses.insert(getAnalysis(n));
    }
    return analyses;
  }


  template <typename Fn>
  inline Fn evil_cast(void* ptr) {
    return reinterpret_cast<Fn>(reinterpret_cast<size_t>(ptr));
  }


  AnalysisBuilders& AnalysisLoader::loadAnalysisBuildersFromFile(const string& filename, AnalysisBuilders& builders) {
    void* handle = dlopen((filename).c_str(), RTLD_LAZY); 
    if (!handle) {
      cerr << "Cannot open " << filename << ": " << dlerror() << endl;
      return builders;
    }
    // Store a list of libs to be neatly closed later
    _handles.insert(handle);

    // Perform dodgy cast to factory function     
    void* fnptr = dlsym(handle, "getAnalysisBuilders");
    anabuilders_fn getBuilders = evil_cast<anabuilders_fn>(fnptr);
    if (!getBuilders) {
      _handles.erase(handle);
      dlclose(handle);
      return builders;
    }
    
    AnalysisBuilders mybuilders = getBuilders();
    foreach (AnalysisBuilders::value_type& b, mybuilders) {
      builders[b.first] = b.second;
    }
    
    return builders;
  }


  void AnalysisLoader::closeAnalysisBuilders() {
    foreach (void* h, _handles) {
      if (h) dlclose(h);
    }
    _handles.clear();
  }
  

  AnalysisBuilders& AnalysisLoader::loadAnalysisBuildersFromDir(const string& dirname, AnalysisBuilders& builders) {
    #ifdef LIB_SUFFIX
    const string libsuffix(LIB_SUFFIX);
    #else
    const string libsuffix(".so");
    #endif

    set<string> libfiles;
    oslink::directory dir(dirname);
    while (dir) {
      string filename = dir.next();
      size_t posn = filename.find(libsuffix);
      if (posn == string::npos || posn != filename.length()-libsuffix.length()) continue;
      if (filename.find("Rivet") == string::npos) continue;
      libfiles.insert(filename);
    }
    
    foreach (const string& l, libfiles) {        
      // Make sure this is an abs path
      /// @todo Sys-dependent path separator instead of "/"
      const string filename  = dirname + "/" + l;
      loadAnalysisBuildersFromFile(filename, builders);        
    }
    
    return builders;
  }
  
  
  AnalysisBuilders& AnalysisLoader::loadAnalysisBuildersFromDirs(const vector<string>& dirnames, AnalysisBuilders& builders) {
    foreach (const string& d, dirnames) {
      loadAnalysisBuildersFromDir(d, builders);
    }
    return builders;
  }
  
  
  void AnalysisLoader::loadAnalyses() {
    vector<string> dirs;
    char* env = 0;
    
    // Always (try to) use the Rivet library install path
    dirs.push_back(getLibPath());
    
    // Then use the Rivet analysis path variable
    env = getenv("RIVET_ANALYSIS_PATH");
    if (env) dirs += split(env);
    
    // And then the user's (non-system) library path
    env = getenv("LD_LIBRARY_PATH");
    if (env) dirs += split(env);
    
    // Use the current dir, and its libtool subdir, too.
    dirs.push_back("."); 
    dirs.push_back(".libs"); 

    
    // Load libs here
    loadAnalysisBuildersFromDirs(dirs, _analysisbuilders);
    _loaded = true;
  }
  

}
