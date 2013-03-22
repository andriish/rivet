#include "Rivet/Particle.hh"
#include "Rivet/RivetBoost.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  bool Particle::hasAncestor(PdgId pdg_id) const {
    /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
    GenVertex* prodVtx = genParticle()->production_vertex();
    if (prodVtx == 0) return false;
    foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
      if (ancestor->pdg_id() == pdg_id) return true;
    }
    return false;
  }


  bool Particle::fromDecay() const {
    /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
    GenVertex* prodVtx = genParticle()->production_vertex();
    if (prodVtx == NULL) return false;
    foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
      const PdgId pid = ancestor->pdg_id();
      if (ancestor->status() == 2 && (PID::isHadron(pid) || abs(pid) == PID::TAU)) return true;
    }
    return false;
  }


}
