// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FinalStateProjection class.
//

#include "FinalStateProjection.h"

using namespace Rivet;

FinalStateProjection::~FinalStateProjection() {}

int FinalStateProjection::cmp(const Projection & p) const {
  return 0;
}

void FinalStateProjection::project(const Event & e) {
  for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
	pi != e.genEvent().particles_end(); ++pi ) {
    if ( (*pi)->status() == 0 ) {
      Particle p;
      p.original = *pi;
      p.id = (*pi)->pdg_id();
      p.momentum = (*pi)->momentum();
      p.mass = (*pi)->generatedMass();
      theParticles.push_back(p);
    }
  }
}

