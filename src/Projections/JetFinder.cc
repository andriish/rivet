// -*- C++ -*-
#include "Rivet/Projections/JetFinder.hh"

namespace Rivet {


  JetFinder::JetFinder(const FinalState& fs, Muons usemuons, Invisibles useinvis)
    : _useMuons(usemuons), _useInvisibles(useinvis)
  {
    setName("JetFinder");
    declare(fs, "FS");
    MSG_TRACE("Constructing final state in JetFinder");
    VisibleFinalState vfs(fs);
    // MSG_DEBUG("Making visible final state from provided FS");
    MSG_TRACE("Declaring visible final state in JetFinder");
    declare(vfs, "VFS");
  }


}
