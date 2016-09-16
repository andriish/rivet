#include "Rivet/Particle.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  Particle& Particle::transformBy(const LorentzTransform& lt) {
    _momentum = lt.transform(_momentum);
    return *this;
  }


  bool Particle::isVisible() const {
    // Charged particles are visible
    if ( PID::threeCharge(pid()) != 0 ) return true;
    // Neutral hadrons are visible
    if ( PID::isHadron(pid()) ) return true;
    // Photons are visible
    if ( pid() == PID::PHOTON ) return true;
    // Gluons are visible (for parton level analyses)
    if ( pid() == PID::GLUON ) return true;
    // Everything else is invisible
    return false;
  }


  bool Particle::isStable() const {
    return genParticle() != NULL &&
      genParticle()->status() == 1 &&
      genParticle()->end_vertex() == NULL;
  }


  vector<Particle> Particle::parents(const Cut& c) const {
    vector<Particle> rtn;
    /// @todo Remove this const mess crap when HepMC doesn't suck
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->production_vertex() );
    if (gv == NULL) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::children))
    //   rtn += Particle(gp);
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::parents); it != gv->particles_end(HepMC::parents); ++it) {
      const Particle p(*it);
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      rtn += p;
    }
    return rtn;
  }


  vector<Particle> Particle::children(const Cut& c) const {
    vector<Particle> rtn;
    if (isStable()) return rtn;
    /// @todo Remove this const mess crap when HepMC doesn't suck
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->end_vertex() );
    if (gv == NULL) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::children))
    //   rtn += Particle(gp);
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::children); it != gv->particles_end(HepMC::children); ++it) {
      const Particle p(*it);
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      rtn += p;
    }
    return rtn;
  }


  /// @todo Insist that the current particle is post-hadronization, otherwise throw an exception?
  /// @todo Use recursion through replica-avoiding functions to avoid bookkeeping duplicates
  vector<Particle> Particle::allDescendants(const Cut& c, bool remove_duplicates) const {
    vector<Particle> rtn;
    if (isStable()) return rtn;
    /// @todo Remove this const mess crap when HepMC doesn't suck
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->end_vertex() );
    if (gv == NULL) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::descendants))
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::descendants); it != gv->particles_end(HepMC::descendants); ++it) {
      const Particle p(*it);
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      if (remove_duplicates && (*it)->end_vertex() != NULL) {
        // size_t n = 0; ///< @todo Only remove 1-to-1 duplicates?
        bool dup = false;
        /// @todo Yuck, HepMC
        for (GenVertex::particle_iterator it2 = (*it)->end_vertex()->particles_begin(HepMC::children); it2 != (*it)->end_vertex()->particles_end(HepMC::children); ++it2) {
          // n += 1; if (n > 1) break;
          if ((*it)->pdg_id() == (*it2)->pdg_id()) { dup = true; break; }
        }
        if (dup) continue;
      }
      rtn += p;
    }
    return rtn;
  }


  /// @todo Insist that the current particle is post-hadronization, otherwise throw an exception?
  vector<Particle> Particle::stableDescendants(const Cut& c) const {
    vector<Particle> rtn;
    if (isStable()) return rtn;
    /// @todo Remove this const mess crap when HepMC doesn't suck
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->end_vertex() );
    if (gv == NULL) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::descendants))
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::descendants); it != gv->particles_end(HepMC::descendants); ++it) {
      // if ((*it)->status() != 1 || (*it)->end_vertex() != NULL) continue;
      const Particle p(*it);
      if (!p.isStable()) continue;
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      rtn += p;
    }
    return rtn;
  }


  double Particle::flightLength() const {
    if (isStable()) return -1;
    if (genParticle() == NULL) return 0;
    if (genParticle()->production_vertex() == NULL) return 0;
    const HepMC::FourVector v1 = genParticle()->production_vertex()->position();
    const HepMC::FourVector v2 = genParticle()->end_vertex()->position();
    return sqrt(sqr(v2.x()-v1.x()) + sqr(v2.y()-v1.y()) + sqr(v2.z()-v1.z()));
  }


  bool Particle::hasParent(PdgId pid) const {
    return _hasRelativeWith(HepMC::parents, hasPID(pid));
  }

  bool Particle::hasParentWith(const Cut& c) const {
    return hasParentWith([&](const Particle& p){return c->accept(p);});
  }


  bool Particle::hasAncestor(PdgId pid) const {
    return _hasRelativeWith(HepMC::ancestors, hasPID(pid));
  }

  bool Particle::hasAncestorWith(const Cut& c) const {
    return hasAncestorWith([&](const Particle& p){return c->accept(p);});
  }


  bool Particle::fromBottom() const {
    return _hasRelativeWith(HepMC::ancestors, [](const Particle& p){
        return p.genParticle()->status() == 2 && p.isHadron() && p.hasBottom();
      });
    // const GenVertexPtr prodVtx = genParticle()->production_vertex();
    // if (prodVtx == NULL) return false;
    // foreach (const GenParticlePtr ancestor, particles(prodVtx, HepMC::ancestors)) {
    //   const PdgId pid = ancestor->pdg_id();
    //   if (ancestor->status() == 2 && (PID::isHadron(pid) && PID::hasBottom(pid))) return true;
    // }
    // return false;
  }


  bool Particle::fromCharm() const {
    return _hasRelativeWith(HepMC::ancestors, [](const Particle& p){
        return p.genParticle()->status() == 2 && p.isHadron() && p.hasCharm();
      });
    // const GenVertexPtr prodVtx = genParticle()->production_vertex();
    // if (prodVtx == NULL) return false;
    // foreach (const GenParticlePtr ancestor, particles(prodVtx, HepMC::ancestors)) {
    //   const PdgId pid = ancestor->pdg_id();
    //   if (ancestor->status() == 2 && (PID::isHadron(pid) && PID::hasCharm(pid) && !PID::hasBottom(pid))) return true;
    // }
    // return false;
  }


  bool Particle::fromHadron() const {
    return _hasRelativeWith(HepMC::ancestors, [](const Particle& p){
        return p.genParticle()->status() == 2 && p.isHadron();
      });
    // const GenVertexPtr prodVtx = genParticle()->production_vertex();
    // if (prodVtx == NULL) return false;
    // foreach (const GenParticlePtr ancestor, particles(prodVtx, HepMC::ancestors)) {
    //   const PdgId pid = ancestor->pdg_id();
    //   if (ancestor->status() == 2 && PID::isHadron(pid)) return true;
    // }
    // return false;
  }


  bool Particle::fromTau(bool prompt_taus_only) const {
    if (prompt_taus_only && fromHadron()) return false;
    return _hasRelativeWith(HepMC::ancestors, [](const Particle& p){
        return p.genParticle()->status() == 2 && isTau(p);
      });
    // const GenVertexPtr prodVtx = genParticle()->production_vertex();
    // if (prodVtx == NULL) return false;
    // foreach (const GenParticlePtr ancestor, particles(prodVtx, HepMC::ancestors)) {
    //   const PdgId pid = ancestor->pdg_id();
    //   if (ancestor->status() == 2 && abs(pid) == PID::TAU) return true;
    // }
    // return false;
  }


  // bool Particle::fromDecay() const {
  //   const GenVertexPtr prodVtx = genParticle()->production_vertex();
  //   if (prodVtx == NULL) return false;
  //   foreach (const GenParticlePtr ancestor, particles(prodVtx, HepMC::ancestors)) {
  //     const PdgId pid = ancestor->pdg_id();
  //     if (ancestor->status() == 2 && (PID::isHadron(pid) || abs(pid) == PID::TAU)) return true;
  //   }
  //   return false;
  // }


  bool Particle::isPrompt(bool from_prompt_tau, bool from_prompt_mu) const {
    if (genParticle() == NULL) return false; // no HepMC connection, give up! Throw UserError exception?
    const GenVertexPtr prodVtx = genParticle()->production_vertex();
    if (prodVtx == NULL) return false; // orphaned particle, has to be assume false
    const pair<GenParticlePtr, GenParticlePtr> beams = prodVtx->parent_event()->beam_particles();

    /// @todo Would be nicer to be able to write this recursively up the chain, exiting as soon as a parton or string/cluster is seen
    foreach (const GenParticlePtr ancestor, Rivet::particles(prodVtx, HepMC::ancestors)) {
      const PdgId pid = ancestor->pdg_id();
      if (ancestor->status() != 2) continue; // no non-standard statuses or beams to be used in decision making
      if (ancestor == beams.first || ancestor == beams.second) continue; // PYTHIA6 uses status 2 for beams, I think... (sigh)
      if (PID::isParton(pid)) continue; // PYTHIA6 also uses status 2 for some partons, I think... (sigh)
      if (PID::isHadron(pid)) return false; // prompt particles can't be from hadron decays
      if (abs(pid) == PID::TAU && abspid() != PID::TAU && !from_prompt_tau) return false; // allow or ban particles from tau decays (permitting tau copies)
      if (abs(pid) == PID::MUON && abspid() != PID::MUON && !from_prompt_mu) return false; // allow or ban particles from muon decays (permitting muon copies)
    }
    return true;
  }


  ///////////////////////
  // From Tools/ParticleUtils.hh -- typically to avoid cyclic includes/refs to Cut definition

  FirstParticleWith::FirstParticleWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  FirstParticleWithout::FirstParticleWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  LastParticleWith::LastParticleWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  LastParticleWithout::LastParticleWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }


  Particles& ifilter_select(Particles& particles, const Cut& c) {
    if (c == Cuts::OPEN) return particles;
    // return ifilter_select(particles, *c);
    return ifilter_select(particles, [&](const Particle& p){return c->accept(p);});
  }


  Particles& ifilter_discard(Particles& particles, const Cut& c) {
    if (c == Cuts::OPEN) { particles.clear(); return particles; }
    // return ifilter_discard(particles, *c);
    return ifilter_discard(particles, [&](const Particle& p){return c->accept(p);});
  }


  ///////////////////////



  string to_str(const Particle& p) {
    string pname;
    try {
      pname = PID::toParticleName(p.pid());
    } catch (...) {
      pname = "PID=" + to_str(p.pid());
    }
    stringstream out;
    out << pname << " @ " << p.momentum() << " GeV";
    return out.str();
  }


  string to_str(const ParticlePair& pair) {
    stringstream out;
    out << "[" << pair.first << ", " << pair.second << "]";
    // out << "["
    //     << PID::toParticleName(pair.first.pid()) << " @ "
    //     << pair.first.momentum().E()/GeV << " GeV, "
    //     << PID::toParticleName(pair.second.pid()) << " @ "
    //     << pair.second.momentum().E()/GeV << " GeV]";
    return out.str();
  }



}
