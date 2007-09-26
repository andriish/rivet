// -*- C++ -*-
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/Utils.hh"
#include "osdir.hh"
#include <dlfcn.h>


namespace Rivet {

  bool AnalysisLoader::_loaded = false;
  AnalysisBuilders AnalysisLoader::_analysisbuilders;
  set<void*> AnalysisLoader::_handles;


  set<string> AnalysisLoader::getAllAnalysisNames() {
    if (!_loaded) loadAnalyses();
    set<string> names;
    for (AnalysisBuilders::const_iterator ab = _analysisbuilders.begin(); 
         ab != _analysisbuilders.end(); ++ab) {
      names.insert(ab->first);
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
    for (set<string>::const_iterator n = names.begin(); n != names.end(); ++n) {
      analyses.insert(getAnalysis(*n));
    }
    return analyses;
  }


  AnalysisBuilders& AnalysisLoader::loadAnalysisBuildersFromFile(const string& filename, AnalysisBuilders& builders) {      
    //cout << "Trying " << filename << endl;
    void* handle = dlopen((filename).c_str(), RTLD_LAZY); 
    if (!handle) {
      cerr << "Cannot open " << filename << ": " << dlerror() << endl;
      return builders;
    }
    // Store a list of libs to be neatly closed later
    _handles.insert(handle);
     
    
    anabuilders_fn getBuilders = (anabuilders_fn) dlsym(handle, "getAnalysisBuilders");
    if (!getBuilders) {
      _handles.erase(handle);
      dlclose(handle);
      return builders;
    }
    
    AnalysisBuilders mybuilders = getBuilders();
    for (AnalysisBuilders::iterator b = mybuilders.begin(); b != mybuilders.end(); ++b) {
      builders[b->first] = b->second;
    }
    
    return builders;
  }


  void AnalysisLoader::closeAnalysisBuilders() {
    for (set<void*>::iterator h = _handles.begin(); h != _handles.end(); ++h) {
      if (*h) dlclose(*h);
    }
    _handles.clear();
  }
    

  AnalysisBuilders& AnalysisLoader::loadAnalysisBuildersFromDir(const string& dirname, AnalysisBuilders& builders) {
    set<string> libfiles;
    oslink::directory dir(dirname);
    while (dir) {
      string filename = dir.next();
      
      // Require that lib filename matches "*Rivet*.{so,dylib,dll}" (basic glob)
      // @todo Use #define SYSDSO ".so", ".dylib", ".dll"
      #define SYSDSO string(".so")
      size_t posn = filename.find(SYSDSO);
      if (posn == string::npos || posn != filename.length()-SYSDSO.length()) continue;
      if (filename.find("Rivet") == string::npos) continue;
      libfiles.insert(filename);
    }
    
    for (set<string>::const_iterator l = libfiles.begin(); l != libfiles.end(); ++l) {        
      // Make sure this is an abs path
      /// @todo Sys-dependent path separator instead of "/"
      const string filename  = dirname + "/" + *l;
      loadAnalysisBuildersFromFile(filename, builders);        
    }
    
    return builders;
  }
  
  
  AnalysisBuilders& AnalysisLoader::loadAnalysisBuildersFromDirs(const vector<string>& dirnames, AnalysisBuilders& builders) {
    for (vector<string>::const_iterator d = dirnames.begin(); d != dirnames.end(); ++d) {
      loadAnalysisBuildersFromDir(*d, builders);
    }
    return builders;
  }
  
  
  void AnalysisLoader::loadAnalyses() {
    vector<string> dirs;
    char* env = 0;
    
    // Always (try to) use the Rivet library install path
    //cout << "** LIBS = " << getInstalledLibPath() << endl;
    dirs.push_back(getInstalledLibPath());
    
    // Then use the Rivet analysis path variable
    env = getenv("RIVET_ANALYSIS_PATH");
    if (env) dirs += split(env);
    
    // And then the user's (non-system) library path
    env = getenv("LD_LIBRARY_PATH");
    if (env) dirs += split(env);
    
    // Use the current dir, too.
    dirs.push_back(".");      
    
    // Load libs here
    loadAnalysisBuildersFromDirs(dirs, _analysisbuilders);
    _loaded = true;
  }
  

}
