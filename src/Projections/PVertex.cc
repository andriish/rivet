// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/PVertex.hh"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"

namespace Rivet {


  void PVertex::project(const Event& e) {
    // We'll *try* to do it right, in case any generators are doing the right thing...
    _thePVertex = e.genEvent().signal_process_vertex();
    getLog() << Log::DEBUG << "PVertex ptr from HepMC = " << _thePVertex << endl;
    if (!_thePVertex) {
      // Since no signal vertices are filled in existing Fortran & C++ MC's,
      // the decay vertex from first particle (supposed to be beam particle)
      // is set as Primary Vertex to get right vertex positions. 
      _thePVertex = (*(e.genEvent().particles_begin()))->end_vertex();
    }
    assert(_thePVertex);
    const unsigned int pVertexParticleSize = _thePVertex->particles_in_size();
    if (pVertexParticleSize != 2 ) {
      throw Error("Wrong number of Primary Vertex particles: " + pVertexParticleSize);
    }
  }


}
