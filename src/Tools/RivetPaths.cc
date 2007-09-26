#include "Rivet/Rivet.hh"

#define STRINGIFY(x) #x

namespace Rivet {

  const string getInstalledDataPath() {
    return RIVETDATADIR;
  }

  const string getInstalledLibPath() {
    return RIVETLIBDIR;
  }

}
