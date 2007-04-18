// -*- C++ -*-

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;


int FinalState::compare(const Projection & p) const {
  const FinalState & other =
    dynamic_cast<const FinalState &>(p);
  return cmp(etamin, other.etamin) || cmp(etamax, other.etamax) ||
    cmp(ptmin, other.ptmin);
}

void FinalState::project(const Event & e) {
  theParticles.clear();
  for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
        pi != e.genEvent().particles_end(); ++pi ) {
    // Only include particles which are final state (status = 1) and which
    // pass the eta and phi cuts. The eta cut is pre-tested by checking if the
    // x and y components of the momentum are non-zero since CLHEP might throw
    // an exception otherwise.
    if ( (*pi)->status() == 1 &&
         fabs((*pi)->momentum().x()) > 0.0 &&
         fabs((*pi)->momentum().y()) > 0.0 &&
         (*pi)->momentum().eta() > etamin &&
         (*pi)->momentum().eta() < etamax &&
         (*pi)->momentum().perp() >= ptmin )
      theParticles.push_back(Particle(**pi));
  }
}

