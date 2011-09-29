// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  void CentralEtHCM::project(const Event& e) {
    const DISFinalState& fs = applyProjection<DISFinalState>(e, "FS");
    _sumet = 0.0;
    foreach (const Particle& p, fs.particles()) {
      /// @todo Generalise rapidity cut value
      if (fabs(p.momentum().rapidity()) < 0.5) _sumet += p.momentum().Et();
    }
  }


}
