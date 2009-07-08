#include "Rivet/Rivet.hh"
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

  const string getRivetgunDataPath() {
    return getDataPath() + "/AGILe";
  }


}
