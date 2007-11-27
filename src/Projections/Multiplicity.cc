// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Cmp.hh"
#include "HepPDT/ParticleID.hh"


namespace Rivet {

  int Multiplicity::compare(const Projection& p) const {
    const Multiplicity& other = dynamic_cast<const Multiplicity&>(p);
    return pcmp(_fsproj, other._fsproj);
  }


  void Multiplicity::project(const Event& e) {
    // Project into final state
    const FinalState& fs = e.applyProjection(_fsproj);

    // Increment counters
    _totalMult = fs.particles().size();
    _hadMult = 0;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      const HepPDT::ParticleID pInfo( p->getPdgId() );
      const bool isHadron = pInfo.isHadron();
      if (isHadron) ++_hadMult;
    }
  }

}
