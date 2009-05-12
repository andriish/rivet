// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  void Beam::project(const Event& e) {
    assert(e.genEvent().particles_size() >= 2);
    std::pair<HepMC::GenParticle*, HepMC::GenParticle*> beams = 
      e.genEvent().beam_particles();
    assert(beams.first);
    _theBeams.first = *(beams.first);
    assert(beams.second);
    _theBeams.second = *(beams.second);

    getLog() << Log::DEBUG << "Beam particle IDs = " 
             << _theBeams.first.pdgId() << ", "
             << _theBeams.second.pdgId() << endl;
  }


  const double Beam::sqrtS() const {
    const double mom1 = beams().first.momentum().pz();
    const double mom2 = beams().second.momentum().pz();
    assert(sign(mom1) != sign(mom2));
    double sqrts = 0.0;
    if (fuzzyEquals(fabs(mom1), fabs(mom2))) {
      sqrts = fabs(mom1) + fabs(mom2);
    } else {
      /// @todo Implement general sqrt(s) for asymmetric beams... requires particle masses.
      throw Error("Asymmetric beams... calculation of sqrt(S) not yet implemented");
    }
    getLog() << Log::DEBUG << "sqrt(s) = " << sqrts/GeV << " GeV" << endl;
    return sqrts;
  }


  /////////////////////////////////////////////////


  ParticlePair beams(const Event& e) {
    Beam beamproj;
    beamproj.project(e);
    return beamproj.beams();
  }

  BeamPair beamIds(const Event& e) {
    Beam beamproj;
    beamproj.project(e);
    return beamproj.beamIDs();
  }


}
