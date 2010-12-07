#include "Rivet/Rivet.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/RivetBoost.hh"
#include "binreloc.h"

namespace Rivet {


  const string getLibPath() {
    BrInitError error;
    br_init_lib(&error);
    char* temp = br_find_lib_dir(DEFAULTLIBDIR);
    const string libdir(temp);
    free (temp);
    return libdir;
  }

  const string getDataPath() {
    BrInitError error;
    br_init_lib(&error);
    char* temp = br_find_data_dir(DEFAULTDATADIR);
    const string sharedir(temp);
    free (temp);
    return sharedir;
  }

  const string getRivetDataPath() {
    return getDataPath() + "/Rivet";
  }


  const vector<string> getAnalysisLibPaths() {
    vector<string> dirs;
    char* env = 0;
    env = getenv("RIVET_ANALYSIS_PATH");
    if (env) {
      // Use the Rivet analysis path variable if set...
      dirs += pathsplit(env);
    } else {
      // ... otherwise fall back to the Rivet library install path
      dirs += getLibPath();
    }
    return dirs;
  }


  void setAnalysisLibPaths(const vector<string>& paths) {
    const string pathstr = pathjoin(paths);
    setenv("RIVET_ANALYSIS_PATH", pathstr.c_str(), 1);
  }


  void addAnalysisLibPath(const string& extrapath) {
    vector<string> paths = getAnalysisLibPaths();
    paths.push_back(extrapath);
    setAnalysisLibPaths(paths);
  }


  const vector<string> getAnalysisRefPaths() {
    vector<string> dirs;
    char* env = 0;
    env = getenv("RIVET_REF_PATH");
    if (env) {
      // Use the Rivet analysis path variable if set...
      dirs += pathsplit(env);
    } else {
      // ... otherwise fall back to the Rivet data install path
      dirs += getRivetDataPath();
      // ... and also add any analysis plugin search dirs for convenience
      dirs += getAnalysisLibPaths();
    }
    return dirs;
  }


  const vector<string> getAnalysisInfoPaths() {
    vector<string> dirs;
    char* env = 0;
    env = getenv("RIVET_INFO_PATH");
    if (env) {
      // Use the Rivet analysis path variable if set...
      dirs += pathsplit(env);
    } else {
      // ... otherwise fall back to the Rivet data install path
      dirs += getRivetDataPath();
      // ... and also add any analysis plugin search dirs for convenience
      dirs += getAnalysisLibPaths();
    }
    return dirs;
  }


}
