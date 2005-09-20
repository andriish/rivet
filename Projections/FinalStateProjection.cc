// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FinalStateProjection class.
//

#include "FinalStateProjection.h"

using namespace Rivet;

FinalStateProjection::~FinalStateProjection() {}

int FinalStateProjection::cmp(const Projection & p) const {
  const FinalStateProjection & other =
    dynamic_cast<const FinalStateProjection &>(p);
  if ( etamin < other.etamin ) return -1;
  else if ( etamin > other.etamin ) return 1;
  if ( etamax < other.etamax ) return -1;
  else if ( etamax > other.etamax ) return 1;
  return 0;
}

void FinalStateProjection::project(const Event & e) {
  for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
	pi != e.genEvent().particles_end(); ++pi ) {
    if ( (*pi)->status() == 1 &&
	 (*pi)->momentum().eta() > etamin &&
	 (*pi)->momentum().eta() < etamax)
      theParticles.push_back(Particle(**pi));
  }
}

