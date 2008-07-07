// -*- C++ -*-
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {

  int FinalState::compare(const Projection& p) const {
    const FinalState& other = dynamic_cast<const FinalState&>(p);
    return \
      cmp(_etamin, other._etamin) || 
      cmp(_etamax, other._etamax) || 
      cmp(_ptmin, other._ptmin);
  }


  void FinalState::project(const Event& e) {
    Log& log = getLog();
    _theParticles.clear();

    for (GenEvent::particle_const_iterator p = e.genEvent().particles_begin();
         p != e.genEvent().particles_end(); ++p) {
      // Only include particles which are final state (status = 1) and which
      // pass the eta and phi cuts. The eta cut is pre-tested by checking if the
      // x and y components of the momentum are non-zero since the vectors might
      // throw an exception otherwise.
      const bool passed = accept(**p);  //GIULIO
      if (log.isActive(Log::TRACE)) {
        log << Log::TRACE << std::boolalpha 
            << "ID = " << (*p)->pdg_id() << ", status = " << (*p)->status() << ", pT = " << (*p)->momentum().perp() 
            << ", eta = " << (*p)->momentum().eta() << ": result = " << passed << endl;
      }
      if (passed) _theParticles.push_back(Particle(**p));
    }
    log << Log::DEBUG << "Number of final-state particles = " 
        << _theParticles.size() << endl;
  }

}
