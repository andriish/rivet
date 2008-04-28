// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  int TotalVisibleMomentum::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }


  void TotalVisibleMomentum::project(const Event& e) {
    _momentum = FourMomentum();
    _momentum.px(0.0).py(0.0).pz(0.0).E(0.0);
    _set = 0.0;

    // Project into final state
    const FinalState& fs = applyProjection<FinalState>(e, "FS");

    // Get hadron and charge info for each particle, and fill counters appropriately
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      _momentum += p->getMomentum();
      _set += pT(p->getMomentum());
    }
  }

}
