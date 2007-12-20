// -*- C++ -*-

#include "Rivet/Projections/WZandh.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"


namespace Rivet {

  int WZandh::compare(const Projection & p) const {
    // They have to be same because there are no parameters.
    return true;
  }

  /// @todo Define Z W and h more carefully/model independently, based on final state.
  void WZandh::project(const Event & e) {
    _theWens.clear();
    _theZees.clear();
    _theWmns.clear();
    _theZmms.clear();
    _theWtns.clear();
    _theZtts.clear();
    _theZnns.clear();
    _theWqqs.clear();
    _theZqqs.clear();
    _thehs.clear();
    for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
          pi != e.genEvent().particles_end(); ++pi ) {
      if ( (*pi)->status() != 1 ) {
        const int id = abs((*pi)->pdg_id());
        /// @todo Use ParticleName enum values for clarity.
        if (id == 23 || id == 24 || id == 25 ){
          // This is a W, Z or h.
          // Now find out whether it is the last one before the decay.
          // @todo Trace the decay products of the bosons properly.
          if (id == 23) {
            _theZees.push_back(Particle(**pi));
          } else if (id == 24) {
            _theWens.push_back(Particle(**pi));	
          } else if (id == 25) {
            _thehs.push_back(Particle(**pi));	
          }         
        }
      }
    }
  }
  
}
