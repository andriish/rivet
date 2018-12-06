#include "Rivet/Jet.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  Jet& Jet::clear() {
    _momentum = FourMomentum();
    _pseudojet.reset(0,0,0,0);
    _particles.clear();
    return *this;
  }


  Jet& Jet::setState(const FourMomentum& mom, const Particles& particles, const Particles& tags) {
    clear();
    _momentum = mom;
    _pseudojet = fastjet::PseudoJet(mom.px(), mom.py(), mom.pz(), mom.E());
    _particles = particles;
    _tags = tags;
    return *this;
  }


  Jet& Jet::setState(const fastjet::PseudoJet& pj, const Particles& particles, const Particles& tags) {
    clear();
    _pseudojet = pj;
    _momentum = FourMomentum(pj.e(), pj.px(), pj.py(), pj.pz());
    _particles = particles;
    _tags = tags;
    // if (_particles.empty()) {
    //   for (const fastjet::PseudoJet pjc : _pseudojet.constituents()) {
    //     // If there is no attached user info, we can't create a meaningful particle, so skip
    //     if (!pjc.has_user_info<RivetFJInfo>()) continue;
    //     const RivetFJInfo& fjinfo = pjc.user_info<RivetFJInfo>();
    //     // Don't add ghosts to the particles list
    //     if (fjinfo.isGhost) continue;
    //     // Otherwise construct a Particle from the PseudoJet, preferably from an associated GenParticle
    //     ?if (fjinfo.genParticle != NULL) {
    //       _particles.push_back(Particle(fjinfo.genParticle));
    //     } else {
    //       if (fjinfo.pid == 0) continue; // skip if there is a null PID entry in the FJ info
    //       const FourMomentum pjcmom(pjc.e(), pjc.px(), pjc.py(), pjc.pz());
    //       _particles.push_back(Particle(fjinfo.pid, pjcmom));
    //     }
    //   }
    // }
    return *this;
  }


  Jet& Jet::setParticles(const Particles& particles) {
    _particles = particles;
    return *this;
  }


  bool Jet::containsParticle(const Particle& particle) const {
    #if HEPMC_VERSION_CODE >= 3000000
    const int idcode = particle.genParticle()->id();
    for (const Particle& p : particles()) {
      if (p.genParticle()->id() == idcode) return true;
    }
    #else
    const int barcode = particle.genParticle()->barcode();
    for (const Particle& p : particles()) {
      if (p.genParticle()->barcode() == barcode) return true;
    }
    #endif
    return false;
  }


  bool Jet::containsParticleId(PdgId pid) const {
    for (const Particle& p : particles()) {
      if (p.pid() == pid) return true;
    }
    return false;
  }


  bool Jet::containsParticleId(const vector<PdgId>& pids) const {
    for (const Particle& p : particles()) {
      for (PdgId pid : pids) {
        if (p.pid() == pid) return true;
      }
    }
    return false;
  }


  /// @todo Jet::containsMatch(Matcher m) { ... if m(pid) return true; ... }


  Jet& Jet::transformBy(const LorentzTransform& lt) {
    _momentum = lt.transform(_momentum);
    for (Particle& p : _particles) p.transformBy(lt);
    for (Particle& t : _tags) t.transformBy(lt);
    _pseudojet.reset(_momentum.px(), _momentum.py(), _momentum.pz(), _momentum.E()); //< lose ClusterSeq etc.
    return *this;
  }


  double Jet::neutralEnergy() const {
    double e_neutral = 0.0;
    for (const Particle& p : particles()) {
      if (p.charge3() == 0) e_neutral += p.E();
    }
    return e_neutral;
  }


  double Jet::hadronicEnergy() const {
    double e_hadr = 0.0;
    for (const Particle& p : particles()) {
      if (isHadron(p)) e_hadr += p.E();
    }
    return e_hadr;
  }


  bool Jet::containsCharm(bool include_decay_products) const {
    for (const Particle& p : particles()) {
      if (p.abspid() == PID::CQUARK) return true;
      if (isHadron(p) && hasCharm(p)) return true;
      if (include_decay_products) {
        ConstGenVertexPtr gv = p.genParticle()->production_vertex();
        if (gv) {
          for (ConstGenParticlePtr pi : Rivet::particles(gv, Relatives::ANCESTORS)) {
            const PdgId pid2 = pi->pdg_id();
            if (PID::isHadron(pid2) && PID::hasCharm(pid2)) return true;
          }
        }
      }
    }
    return false;
  }


  bool Jet::containsBottom(bool include_decay_products) const {
    for (const Particle& p : particles()) {
      if (p.abspid() == PID::BQUARK) return true;
      if (isHadron(p) && hasBottom(p)) return true;
      if (include_decay_products) {
        ConstGenVertexPtr gv = p.genParticle()->production_vertex();
        if (gv) {
          for (ConstGenParticlePtr pi : Rivet::particles(gv, Relatives::ANCESTORS)) {
            const PdgId pid2 = pi->pdg_id();
            if (PID::isHadron(pid2) && PID::hasBottom(pid2)) return true;
          }
        }
      }
    }
    return false;
  }


  Particles Jet::tags(const Cut& c) const {
    return filter_select(tags(), c);
  }

  Particles Jet::bTags(const Cut& c) const {
    Particles rtn;
    for (const Particle& tp : tags()) {
      if (hasBottom(tp) && c->accept(tp)) rtn.push_back(tp);
    }
    return rtn;
  }

  Particles Jet::cTags(const Cut& c) const {
    Particles rtn;
    for (const Particle& tp : tags()) {
      /// @todo Is making b and c tags exclusive the right thing to do?
      if (hasCharm(tp) && !hasBottom(tp) && c->accept(tp)) rtn.push_back(tp);
    }
    return rtn;
  }

  Particles Jet::tauTags(const Cut& c) const {
    Particles rtn;
    for (const Particle& tp : tags()) {
      if (isTau(tp) && c->accept(tp)) rtn.push_back(tp);
    }
    return rtn;
  }


}
