// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISLepton class.
//

#include "DISLepton.h"

using namespace Rivet;

DISLepton::~DISLepton() {}


void DISLepton::project(const Event & e) {
  vector<GenParticle*> inc =
    e.genEvent().signal_process_vertex()->listParents();
  for ( int i = 0, N = inc.size(); i < N; ++i )
    if ( inc[i]->pdg_id() == idin ) incoming = Particle(*inc[i]);
  double emax = 0.0;
  for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
	pi != e.genEvent().particles_end(); ++pi ) {
    if ( (*pi)->momentum().e() > emax ) {
      // *** ATTENTION *** This is probably not the way to select the
      // scattered lepton
      emax = (*pi)->momentum().e();
      outgoing = Particle(**pi);
    }
  }
}

int DISLepton::cmp(const Projection & p) const {
  const DISLepton & other = dynamic_cast<const DISLepton &>(p);
  if ( idin < other.idin ) return -1;
  else if ( idin > other.idin ) return 1;
  if ( idout < other.idout ) return -1;
  else if ( idout > other.idout ) return 1;
  return 0;
}

