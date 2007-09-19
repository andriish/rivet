// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Cmp.hh"
#include "HepPDT/ParticleID.hh"

namespace Rivet {

  int ChargedFinalState::compare(const Projection& p) const {
    const ChargedFinalState& other = dynamic_cast<const ChargedFinalState&>(p);
    return pcmp(_fsproj, other._fsproj);
  }
  
  
  void ChargedFinalState::project(const Event& e) {
    Log log = getLog();
    const FinalState& fs = e.applyProjection(_fsproj);
    _theParticles.clear();
    const ParticleVector& fsps = fs.particles();
    _theParticles.reserve(fsps.size());
    for (ParticleVector::const_iterator p = fsps.begin(); p != fsps.end(); ++p) {
      HepPDT::ParticleID pInfo = p->getPdgId();
      bool isHadron = pInfo.isHadron();
      if (pInfo.threeCharge() != 0) {
        _theParticles.push_back(*p);
      }
    }
  } 
  
}
