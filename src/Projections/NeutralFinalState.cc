// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Cmp.hh"
#include <algorithm>

namespace Rivet {


  int NeutralFinalState::compare(const Projection& p) const {
    /// @todo: This needs to be fixed!! We can't just compare the FinalStates, we also have to check the Etmin!!
    return mkNamedPCmp(p, "FS");
  }


  void NeutralFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    foreach (const Particle& p, fs.particles()){
      if ((PID::threeCharge(p.pdgId()) != 0) && (p.momentum().Et() > _Etmin)) {
        _theParticles.push_back(p);
        if (getLog().isActive(Log::TRACE)) {
          getLog() << Log::TRACE
                   << "Selected: ID = " << p.pdgId()
                   << ", Et = " << p.momentum().Et()
                   << ", eta = " << p.momentum().eta()
                   << ", charge = " << PID::threeCharge(p.pdgId())/3.0 << endl;
        }
      }
    }
    getLog() << Log::DEBUG << "Number of neutral final-state particles = "
             << _theParticles.size() << endl;
  }


}
