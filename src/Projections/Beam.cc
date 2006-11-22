// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Beam class.
//

#include "Rivet/Projections/Beam.hh"
#include "HepMC/GenVertex.h"
#include <stdexcept>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Beam.tcc"
#endif

using namespace Rivet;
using std::runtime_error;

Beam::~Beam() {}

void Beam::project(const Event& e) {
  //vector<GenParticle*> inc = e.genEvent().signal_process_vertex()->listParents();
  GenVertex* sigvertex = e.genEvent().signal_process_vertex();
  if ( sigvertex->particles_in_size() != 2 ) {
    throw std::runtime_error("Wrong number of beams.");
    /// @todo Maybe we should have our own exception classes?
  }
  // Why can't HepMC just give us the list of particles? *sigh*
  GenVertex::particles_in_const_iterator pp = sigvertex->particles_in_const_begin();
  theBeams.first = Particle(**pp);
  ++pp;
  theBeams.second = Particle(**pp);
}

int Beam::compare(const Projection &) const {
  return 0;
}
