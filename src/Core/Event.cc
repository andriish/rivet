// -*- C++ -*-
#include "Rivet/Event.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {

  double Event::weight() const {
    return genEvent()->weights().empty() ? 1.0 : _genevent.weights()[0];
  }

  ParticlePair Event::beams() const { return Rivet::beams(*this); }

  double Event::sqrtS() const { return Rivet::sqrtS(beams()); }

  double Event::asqrtS() const { return Rivet::asqrtS(beams()); }



  void Event::_init(const GenEvent& ge) {
    // Use Rivet's preferred units if possible
    #ifdef HEPMC_HAS_UNITS
    _genevent.use_units(HepMC::Units::GEV, HepMC::Units::MM);
    #endif
  }


  const Particles& Event::allParticles() const {
    if (_particles.empty()) { //< assume that empty means no attempt yet made
      for (ConstGenParticlePtr gp : HepMCUtils::particles(genEvent())) {
        _particles += Particle(gp);
      }
    }
    return _particles;
  }


}
