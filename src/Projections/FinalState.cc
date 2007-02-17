// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FinalState class.
//

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;

FinalState::~FinalState() {}

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
    if ( (*pi)->status() == 1 &&
         (*pi)->momentum().eta() > etamin &&
         (*pi)->momentum().eta() < etamax &&
         (*pi)->momentum().perp() >= ptmin )
      theParticles.push_back(Particle(**pi));
  }
}

