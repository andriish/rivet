// -*- C++ -*-

#include "Rivet/Projections/Beam.hh"
#include "HepMC/GenVertex.h"
#include <stdexcept>

using namespace Rivet;
using std::runtime_error;

void Beam::project(const Event& e) {
  //vector<GenParticle*> inc = e.genEvent().signal_process_vertex()->listParents();
  GenVertex* sigvertex = e.genEvent().signal_process_vertex();
  if ( sigvertex->particles_in_size() != 2 ) {
    throw std::runtime_error("Wrong number of beams.");
    /// @todo Maybe we should have our own exception classes?
  }
  // Why can't HepMC just give us the list of particles? *sigh*
  //ls: Oh it can (see bottom)


  GenVertex::particles_in_const_iterator pp = sigvertex->particles_in_const_begin();
  theBeams.first = Particle(**pp);
  ++pp;
  theBeams.second = Particle(**pp);

  /*
  const HepMC::GenEvent ge = e.genEvent();

  GenEvent::particle_const_iterator pp = ge.particles_begin();
  theBeams.first = Particle(**pp);
  //cout << "first Beam particle id=" << (**pp).pdg_id() << endl;
  ++pp;
  theBeams.second = Particle(**pp);
  //cout << "second Beam particle id=" << (**pp).pdg_id() << endl;
  */


}


int Beam::compare(const Projection &) const {
  return 0;
}
