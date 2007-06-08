// -*- C++ -*-

#include "Rivet/Projections/PVertex.hh"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"
#include <stdexcept>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Beam.tcc"
#endif

using namespace Rivet;
using std::runtime_error;

void PVertex::project(const Event& e) {

  thePVertex = e.genEvent().signal_process_vertex();
  if ( !thePVertex || thePVertex->particles_in_size() != 2 )
    throw std::runtime_error("Wrong number of Primary Vertex particles.");

}

int PVertex::compare(const Projection &) const {
  return 0;
}
