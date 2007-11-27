#include "Rivet/Rivet.hh"
#include "Rivet/Tools/binreloc.h"

namespace Rivet {

  const string getInstalledDataPath() {
    BrInitError error;
    br_init_lib(&error);
    const string sharedir = br_find_data_dir(RIVETDATADIR);
    return sharedir;
  }
  
  const string getInstalledLibPath() {
    BrInitError error;
    br_init_lib(&error);
    const string libdir = br_find_data_dir(RIVETLIBDIR);
    return libdir;
  }

}
