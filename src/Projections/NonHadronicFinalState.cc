// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/NonHadronicFinalState.hh"
#include "Rivet/Cmp.hh"
#include <algorithm>

namespace Rivet {


  int NonHadronicFinalState::compare(const Projection& p) const {
    return FinalState::compare(p);
  }


  bool nonHadronFilter(const Particle& p) {
    return PID::isHadron(p.pdgId());
  }

  void NonHadronicFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(),
                        std::back_inserter(_theParticles), nonHadronFilter);
    MSG_DEBUG("Number of hadronic final-state particles = "
             << _theParticles.size());
  }

}
