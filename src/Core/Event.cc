// -*- C++ -*-
#include "Rivet/Event.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Projections/Beam.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {

  double Event::weight() const {
    static ssize_t WEIGHT_INDEX = -1;
    if (WEIGHT_INDEX < 0) {
      WEIGHT_INDEX = getEnvParam<size_t>("RIVET_WEIGHT_INDEX", 0);
      Log::getLog("Core.Weight") << Log::TRACE << "Got weight index from env/default = "<< WEIGHT_INDEX << endl;
    }
    return genEvent()->weights().empty() ? 1.0 : _genevent.weights()[WEIGHT_INDEX];
  }

  double Event::centrality() const {
    /// @todo Use direct "centrality" property if using HepMC3
    return genEvent()->heavy_ion() ? genEvent()->heavy_ion()->impact_parameter() : -1;
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
      for (const GenParticle* gp : particles(genEvent())) {
        _particles += Particle(gp);
      }
    }
    return _particles;
  }


}
