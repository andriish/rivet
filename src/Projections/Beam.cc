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


  const double Beam::getSqrtS() const {
    const double mom1 = getBeams().first.getMomentum().pz();
    const double mom2 = getBeams().second.getMomentum().pz();
    assert(sign(mom1) != sign(mom2));
    double sqrts = 0.0;
    if (fuzzyEquals(fabs(mom1), fabs(mom2))) {
      sqrts = fabs(mom1) + fabs(mom2);
    } else {
      /// @todo Implement general sqrt(s) for asymmetric beams.
      throw Error("Asymmetric beams... calculation of sqrt(S) not yet implemented");
    }
    getLog() << Log::DEBUG << "sqrt(s) = " << sqrts/GeV << " GeV" << endl;
    return sqrts;
  }


}
