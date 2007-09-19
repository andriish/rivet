// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "HepMC/GenVertex.h"


namespace Rivet {

  void Beam::project(const Event& e) {
    Log& log = getLog();

    // Assume that the first two particles in the event are the beam particles
    /// @todo This is not a robust assumption! HepMC should support directly accessing the beams.
    HepMC::GenEvent::particle_const_iterator bp = e.genEvent().particles_begin();
    _theBeams.first = **bp;
    ++bp;
    _theBeams.second = **bp;

    log << Log::DEBUG << "Beam particle IDs = " 
        << _theBeams.first.getPdgId() << ", "
        << _theBeams.second.getPdgId() << endl;
  }

}
