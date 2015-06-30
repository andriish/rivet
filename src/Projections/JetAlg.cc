// -*- C++ -*-
#include "Rivet/Projections/JetAlg.hh"

namespace Rivet {


  JetAlg::JetAlg(const FinalState& fs, InvisiblesStrategy useinvis)
    : _useInvisibles(useinvis)
  {
    setName("JetAlg");
    addProjection(fs, "FS");
    VisibleFinalState vfs(fs);
    MSG_DEBUG("Making visible final state from provided FS");
    addProjection(vfs, "VFS");
  }


}
