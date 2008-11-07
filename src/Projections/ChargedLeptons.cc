// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

namespace Rivet {


  int ChargedLeptons::compare(const Projection& other) const {
    return mkNamedPCmp(other, "FS");
  }


  void ChargedLeptons::project(const Event& e) {
    // Reset result
    _theChargedLeptons.clear();

    // Get hadron and charge info for each particle, and fill counters appropriately
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    foreach (const Particle& p, fs.particles()) {
      if (PID::isLepton(p.pdgId())) {
        if (PID::threeCharge(p.pdgId()) != 0) {
          // Put it into the C.L. vector
          _theChargedLeptons.push_back(Particle(p));
        }
      }
    }
    getLog() << Log::DEBUG << "Done" << endl;
  }


}
