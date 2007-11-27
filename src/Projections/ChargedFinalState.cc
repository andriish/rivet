// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Cmp.hh"
#include "HepPDT/ParticleID.hh"
#include <algorithm>


namespace Rivet {

  int ChargedFinalState::compare(const Projection& p) const {
    const ChargedFinalState& other = dynamic_cast<const ChargedFinalState&>(p);
    return FinalState::compare(other);
  }
  
  bool chargedParticleFilter(const Particle& p) {
    const HepPDT::ParticleID pInfo( p.getPdgId() );
    return pInfo.threeCharge() != 0;
  }
  
  void ChargedFinalState::project(const Event& e) {
    Log log = getLog();
    FinalState fsp = static_cast<FinalState>(*this);
    const FinalState& fs = e.applyProjection(fsp);
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(), 
                        std::back_inserter(_theParticles), chargedParticleFilter);
    getLog() << Log::DEBUG << "Number of charged final-state particles = " 
             << _theParticles.size() << endl;
  } 
  
}
