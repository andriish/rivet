// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/LossyFinalState.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include <algorithm>


namespace Rivet {

  int LossyFinalState::compare(const Projection& p) const {
    const LossyFinalState& other = dynamic_cast<const LossyFinalState&>(p);
    const int fscmp = FinalState::compare(other);
    if (fscmp) return fscmp;
    return cmp(_lossFraction, other._lossFraction);
  }
  
  void LossyFinalState::project(const Event& e) {
    Log log = getLog();
    FinalState fsp = static_cast<FinalState>(*this);
    const FinalState& fs = e.applyProjection(fsp);
    getLog() << Log::DEBUG << "Pre-loss number of FS particles = " << fs.particles().size() << endl;
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(), 
                        std::back_inserter(_theParticles), RandomFilter(_lossFraction));
    getLog() << Log::DEBUG << "Filtered number of FS particles = " << _theParticles.size() 
             << " (should be ~" << (1-_lossFraction)*100 << "%)" << endl;
  }
  
}
