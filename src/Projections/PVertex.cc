// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/PVertex.hh"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"


namespace Rivet {

  void PVertex::project(const Event& e) {
    _thePVertex = e.genEvent().signal_process_vertex();
    const size_t pVertexParticleSize = _thePVertex->particles_in_size();
    if ( !_thePVertex || pVertexParticleSize != 2 )
      throw runtime_error("Wrong number of Primary Vertex particles: " + pVertexParticleSize);
  }

}
