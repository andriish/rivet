// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Cmp.hh"
#include <algorithm>


namespace Rivet {

  int ChargedFinalState::compare(const Projection& p) const {
    const ChargedFinalState& other = dynamic_cast<const ChargedFinalState&>(p);
    //return FinalState::compare(other);
    return pcmp(*_fsproj, *other._fsproj);
  }

  
  bool chargedParticleFilter(const Particle& p) {
    return PID::threeCharge(p.getPdgId()) == 0;
  }

  
  void ChargedFinalState::project(const Event& e) {
    Log log = getLog();
    /// @todo This goes out of scope immediately!
    //FinalState fsp = static_cast<FinalState>(*this);
    //const FinalState& fs = e.applyProjection(fsp);
    const FinalState& fs = e.applyProjection(*_fsproj);
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(), 
                        std::back_inserter(_theParticles), chargedParticleFilter);
    getLog() << Log::DEBUG << "Number of charged final-state particles = " 
             << _theParticles.size() << endl;
    if (getLog().isActive(Log::TRACE)) {
      for (vector<Particle>::iterator p = _theParticles.begin(); p != _theParticles.end(); ++p) {
        getLog() << Log::TRACE << "Selected: " << p->getPdgId() 
                 << ", charge = " << PID::threeCharge(p->getPdgId())/3.0 << endl;
      }
    }
  } 
  
}
