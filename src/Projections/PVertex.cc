// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/PVertex.hh"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"


namespace Rivet {

  void PVertex::project(const Event& e) {

    _thePVertex = e.genEvent().signal_process_vertex();

    /// since no signal vertices are filled (nor Fortran, neither C++ MC's)
    /// the decay vertex from first particle (supposed to be beam particle)
    /// is set as Primary Vertex to get right vertex positions. 

    //if (_thePVertex) cout << "_thePVertex=0" << endl;
    //_thePVertex = (*(e.genEvent().particles_begin()))->production_vertex();
    _thePVertex = (*(e.genEvent().particles_begin()))->end_vertex();


    const size_t pVertexParticleSize = _thePVertex->particles_in_size();
    //if ( !_thePVertex || pVertexParticleSize != 2 )
    if ( !_thePVertex )
    if (pVertexParticleSize != 2 )
      throw runtime_error("Wrong number of Primary Vertex particles: " + pVertexParticleSize);
  }

}
