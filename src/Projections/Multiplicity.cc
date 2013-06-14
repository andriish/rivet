// -*- C++ -*-
#include "Rivet/Projections/Multiplicity.hh"

namespace Rivet {


  int Multiplicity::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }


  void Multiplicity::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _totalMult = fs.particles().size();
    _hadMult = 0;
    for (Particles::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      if (PID::isHadron(p->pdgId())) ++_hadMult;
    }
  }

}
