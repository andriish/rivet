// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/LossyFinalState.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include <algorithm>

namespace Rivet {


  int LossyFinalState::compare(const Projection& p) const {
    const LossyFinalState& other = pcast<LossyFinalState>(p);
    const int fscmp = mkNamedPCmp(other, "FS");
    if (fscmp) return fscmp;
    return cmp(_lossFraction, other._lossFraction);
  }
  

  void LossyFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    getLog() << Log::DEBUG << "Pre-loss number of FS particles = " << fs.particles().size() << endl;
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(), 
                        std::back_inserter(_theParticles), RandomFilter(_lossFraction));
    getLog() << Log::DEBUG << "Filtered number of FS particles = " << _theParticles.size() 
             << " (should be ~" << (1-_lossFraction)*100 << "%)" << endl;
  }

  
}
