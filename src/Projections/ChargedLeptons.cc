// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {

  int ChargedLeptons::compare(const Projection& p) const {
    const ChargedLeptons& other = dynamic_cast<const ChargedLeptons &>(p);
    return pcmp(*_fsproj, *other._fsproj);
  }


  void ChargedLeptons::project(const Event& e) {
    Log& log = getLog();

    // Reset result
    _theChargedLeptons.clear();

    // Project into final state
    const FinalState& fs = e.applyProjection(*_fsproj);

    // Get hadron and charge info for each particle, and fill counters appropriately
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      if (PID::isLepton(p->getPdgId())) {
        if (PID::threeCharge(p->getPdgId()) != 0) {
          // Put it into the C.L. vector
          _theChargedLeptons.push_back(Particle(*p));
        }
      }
    }
    log << Log::DEBUG << "Done" << endl;
  }

}
