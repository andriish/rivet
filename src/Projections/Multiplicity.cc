// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  int Multiplicity::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }


  void Multiplicity::project(const Event& e) {
    // Project into final state
    const FinalState& fs = applyProjection<FinalState>(e, "FS");

    // Increment counters
    _totalMult = fs.particles().size();
    _hadMult = 0;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      if (PID::isHadron(p->getPdgId())) ++_hadMult;
    }
  }

}
