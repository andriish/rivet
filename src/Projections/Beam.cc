// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "HepMC/GenVertex.h"


namespace Rivet {

  void Beam::project(const Event& e) {
    Log& log = getLog();

    GenVertex* sigvertex = e.genEvent().signal_process_vertex();
    if ( sigvertex->particles_in_size() != 2 ) {
      throw std::runtime_error("Wrong number of beams.");
    }

    /// @todo What's the point of this if it's being overwritten later?
    GenVertex::particles_in_const_iterator pp = sigvertex->particles_in_const_begin();
    _theBeams.first = Particle(**pp);
    ++pp;
    _theBeams.second = Particle(**pp);

    /*
    const HepMC::GenEvent& ge = e.genEvent();
    GenEvent::particle_const_iterator pp = ge.particles_begin();
    theBeams.first = Particle(**pp);
    //cout << "first Beam particle id=" << (**pp).pdg_id() << endl;
    ++pp;
    theBeams.second = Particle(**pp);
    //cout << "second Beam particle id=" << (**pp).pdg_id() << endl;
    */

    HepMC::GenEvent::particle_const_iterator bp = e.genEvent().particles_begin();
    _theBeams.first = **bp;
    ++bp;
    _theBeams.second = **bp;
    log << Log::DEBUG << "Beam particle IDs = " 
        << _theBeams.first.getPdgId() << ", "
        << _theBeams.second.getPdgId() << endl;
  }

}
