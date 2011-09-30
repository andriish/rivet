// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Cmp.hh"
#include <algorithm>

namespace Rivet {


  int VisibleFinalState::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }


  // Since we remove inivisibles from the FinalState in project(),
  // we need a filter where invisible --> true
  bool isInvisibleFilter(const Particle& p) {
    // charged particles are visible
    if ( PID::threeCharge( p.pdgId() ) != 0 ) 
      return false;

    // neutral hadrons are visible
    if ( PID::isHadron( p.pdgId() ) ) 
      return false;

    // photons are visible
    if ( p.pdgId() == PHOTON ) 
      return false;

    // gluons are visible (for parton level analyses)
    if ( p.pdgId() == GLUON ) 
      return false;

    // everything else is invisible
    return true;
  }


  void VisibleFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(),
                        std::back_inserter(_theParticles), isInvisibleFilter);
    MSG_DEBUG("Number of visible final-state particles = "
             << _theParticles.size());
  }

}
