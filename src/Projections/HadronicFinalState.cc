// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Cmp.hh"
#include <algorithm>


namespace Rivet {

  int HadronicFinalState::compare(const Projection& p) const {
    const HadronicFinalState& other = dynamic_cast<const HadronicFinalState&>(p);
    return FinalState::compare(other);
  }
  
  bool hadronFilter(const Particle& p) {
    return ! PID::isHadron(p.getPdgId());
  }
  
  void HadronicFinalState::project(const Event& e) {
    Log log = getLog();
    FinalState fsp = static_cast<FinalState>(*this);
    const FinalState& fs = e.applyProjection(fsp);
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(), 
                        std::back_inserter(_theParticles), hadronFilter);
    getLog() << Log::DEBUG << "Number of hadronic final-state particles = " 
             << _theParticles.size() << endl;
  } 
  
}
