// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FinalStateProjection class.
//

#include "Rivet/Projections/FinalStateProjection.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;

FinalStateProjection::~FinalStateProjection() {}

int FinalStateProjection::compare(const Projection & p) const {
  const FinalStateProjection & other =
    dynamic_cast<const FinalStateProjection &>(p);
  return cmp(etamin, other.etamin) || cmp(etamax, other.etamax);
}

void FinalStateProjection::project(const Event & e) {
  theParticles.clear();
  for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
        pi != e.genEvent().particles_end(); ++pi ) {
    if ( (*pi)->status() == 1 &&
         (*pi)->momentum().eta() > etamin &&
         (*pi)->momentum().eta() < etamax)
      theParticles.push_back(Particle(**pi));
  }
}

