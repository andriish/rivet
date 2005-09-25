// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISLepton class.
//

#include "DISLepton.h"
#include "Rivet/Projections/Cmp.h"

using namespace Rivet;

DISLepton::~DISLepton() {}


void DISLepton::project(const Event & e) {
  const PPair & inc = e(beams)();
  if ( inc.first.id == idin ) incoming = inc.first;
  else if ( inc.second.id == idin ) incoming = inc.second;
  else
    throw runtime_error("DISLepton projector could not find the correct beam.");
  // *** ATTENTION *** Should we have our own exception classes?

  double emax = 0.0;
  for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
	pi != e.genEvent().particles_end(); ++pi ) {
    if ( (*pi)->momentum().e() > emax ) {
      // *** ATTENTION *** This is probably not the correct way to
      // select the scattered lepton
      emax = (*pi)->momentum().e();
      outgoing = Particle(**pi);
    }
  }
}

int DISLepton::compare(const Projection & p) const {
  const DISLepton & other = dynamic_cast<const DISLepton &>(p);
  return pcmp(beams, other.beams) ||
    cmp(idin, other.idin) || cmp(idout, other.idout);
}

