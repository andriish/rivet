// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"


namespace Rivet {

  void Beam::project(const Event& e) {
    Log& log = getLog();

    assert(e.genEvent().particles_size() >= 2);
    std::pair<HepMC::GenParticle*, HepMC::GenParticle*> beams = e.genEvent().beam_particles();
    _theBeams.first = *(beams.first);
    _theBeams.second = *(beams.second);

    log << Log::DEBUG << "Beam particle IDs = " 
        << _theBeams.first.getPdgId() << ", "
        << _theBeams.second.getPdgId() << endl;
  }

}
