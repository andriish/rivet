// -*- C++ -*-
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  int FinalState::compare(const Projection& p) const {
    const FinalState & other = dynamic_cast<const FinalState &>(p);
    return cmp(_etamin, other._etamin) || cmp(_etamax, other._etamax) || cmp(_ptmin, other._ptmin);
  }


  void FinalState::project(const Event& e) {
    _theParticles.clear();


    for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
          pi != e.genEvent().particles_end(); ++pi ) {
      // Only include particles which are final state (status = 1) and which
      // pass the eta and phi cuts. The eta cut is pre-tested by checking if the
      // x and y components of the momentum are non-zero since CLHEP might throw
      // an exception otherwise.
      if ( (*pi)->status() == 1 &&
           fabs((*pi)->momentum().x()) > 0.0 &&
           fabs((*pi)->momentum().y()) > 0.0 &&
           (*pi)->momentum().eta() > _etamin &&
           (*pi)->momentum().eta() < _etamax &&
           (*pi)->momentum().perp() >= _ptmin )
        _theParticles.push_back(Particle(**pi));
    }
    getLog() << Log::DEBUG << "Number of final-state particles = " 
             << _theParticles.size() << endl;
  }

}
