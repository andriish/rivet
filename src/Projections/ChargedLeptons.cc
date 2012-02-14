// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  int ChargedLeptons::compare(const Projection& other) const {
    return mkNamedPCmp(other, "ChFS");
  }


  void ChargedLeptons::project(const Event& evt) {
    // Reset result
    _theParticles.clear();

    // Loop over charged particles and fill vector with leptons
    const FinalState& fs = applyProjection<FinalState>(evt, "ChFS");
    foreach (const Particle& p, fs.particles()) {
      if (PID::isLepton(p.pdgId())) {
        _theParticles += Particle(p);
      }
    }
    std::sort(_theParticles.begin(), _theParticles.end(), cmpParticleByPt);
  }


}
