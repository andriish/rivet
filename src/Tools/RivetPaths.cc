#include "Rivet/Rivet.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/RivetBoost.hh"
#include "binreloc.h"

namespace Rivet {


  const string getLibPath() {
    BrInitError error;
    br_init_lib(&error);
    const string libdir = br_find_lib_dir(DEFAULTLIBDIR);
    return libdir;
  }

  const string getDataPath() {
    BrInitError error;
    br_init_lib(&error);
    const string sharedir = br_find_data_dir(DEFAULTDATADIR);
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
      dirs += split(env);
    } else {
      // ... otherwise fall back to the Rivet library install path
      dirs += getLibPath();
    }
    return dirs;
  }


  const vector<string> getAnalysisRefPaths() {
    vector<string> dirs;
    char* env = 0;
    env = getenv("RIVET_REF_PATH");
    if (env) {
      // Use the Rivet analysis path variable if set...
      dirs += split(env);
    } else {
      // ... otherwise fall back to the Rivet data install path
      dirs += getRivetDataPath();
    }
    return dirs;
  }


  const vector<string> getAnalysisInfoPaths() {
    vector<string> dirs;
    char* env = 0;
    env = getenv("RIVET_INFO_PATH");
    if (env) {
      // Use the Rivet analysis path variable if set...
      dirs += split(env);
    } else {
      // ... otherwise fall back to the Rivet data install path
      dirs += getRivetDataPath();
    }
    return dirs;
  }


}
