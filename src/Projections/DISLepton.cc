// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISLepton class.
//

#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;


void DISLepton::project(const Event & e) {
  const ParticlePair & inc = e.applyProjection(*beams)();
  if ( ( idin == -idin && abs(inc.first.id) == abs(idin) ) || inc.first.id ) {
    incoming = inc.first;
  } else {
    throw runtime_error("DISLepton projector could not find the correct beam.");
    /// @todo Should we have our own exception classes?
  }

  double emax = 0.0;
  for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
        pi != e.genEvent().particles_end(); ++pi ) {
    if ( ( idin == -idin && abs((*pi)->pdg_id()) == abs(idout) ) ||
         (*pi)->pdg_id() == idout && (*pi)->momentum().e() > emax ) {
      /// @todo This is probably not the correct way to select the scattered lepton
      emax = (*pi)->momentum().e();
      outgoing = Particle(**pi);
    }
  }
}

int DISLepton::compare(const Projection & p) const {
  const DISLepton & other = dynamic_cast<const DISLepton &>(p);
  return pcmp(*beams, *other.beams) ||
    cmp(idin, other.idin) || cmp(idout, other.idout);
}

RivetInfo DISLepton::getInfo() const {
  return Projection::getInfo() + beams->getInfo();
}

