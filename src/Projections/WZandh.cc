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
    bool e_decay = false;
    bool mu_decay = false;
    bool tau_decay = false;
    bool nu_decay = false;
    bool q_decay = false;
    for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
          pi != e.genEvent().particles_end(); ++pi ) {
      if ( (*pi)->status() != 1 ) {
        const int id = abs((*pi)->pdg_id());
	/// @todo Use ParticleName enum values for clarity.
        if (id == 23 || id == 24 || id == 25 ){
	  //if (id == ZBOSON || id == WPLUSBOSON || id == HIGGS) {
          // This is a W, Z or h.
          // Now find out whether it is the last one before the decay.
          // Trace the decay products of the bosons properly.
	  bool bosondecay = true;
	  for (GenVertex::particles_out_const_iterator dpi 
		 = (*pi)->end_vertex()->particles_out_const_begin(); 
	       dpi != (*pi)->end_vertex()->particles_out_const_end(); ++dpi) {  
	    int did = abs((*dpi)->pdg_id());

	    if (did == id) bosondecay = false; 
	    else {
	      if (did == ELECTRON) e_decay = true;
	      else if (did == MUON) mu_decay = true;
	      else if (did == TAU) tau_decay = true;
	      else if (did == NU_E || did == NU_MU || did == NU_TAU) nu_decay = true;
	      else if (did >= 1 && did <=5) q_decay = true;
	    }
	  }

	  if (!bosondecay) continue;

	  /// @todo Use ParticleName enum values for clarity.
          //if (id == ZBOSON) {
	  if (id == 23) {
	    if (e_decay) _theZees.push_back(Particle(**pi));
	    else if (mu_decay) _theZmms.push_back(Particle(**pi));
	    else if (tau_decay) _theZtts.push_back(Particle(**pi));
	    else if (nu_decay) _theZnns.push_back(Particle(**pi));
	    else if (q_decay) _theZqqs.push_back(Particle(**pi));
	  /// @todo Use ParticleName enum values for clarity.
	  //} else if (id == WPLUSBOSON) {
	  } else if (id == 24) {
            if (e_decay) _theWens.push_back(Particle(**pi));	
            if (mu_decay) _theWmns.push_back(Particle(**pi));	
            if (tau_decay) _theWtns.push_back(Particle(**pi));	
            if (q_decay) _theWqqs.push_back(Particle(**pi));	
	    /// @todo Use ParticleName enum values for clarity.
	    //} else if (id == HIGGS) {
	  } else if (id == 25) {
            _thehs.push_back(Particle(**pi));	
          }         
        }
      }
    }
  }
  
}
