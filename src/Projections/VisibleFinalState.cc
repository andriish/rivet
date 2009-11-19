// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/Utils.hh"
#include <algorithm>

namespace Rivet {


  int VisibleFinalState::compare(const Projection& p) const {
    return mkNamedPCmp(p, "VFS");
  }


  void VisibleFinalState::project(const Event& e) {
    const FinalState& vfs = applyProjection<FinalState>(e, "VFS");
    _theParticles = vfs.particles();
  }


}
