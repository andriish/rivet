// -*- C++ -*-

#include "Rivet/Projections/WZandh.hh"
#include "Rivet/Projections/Cmp.hh"
#include "HepPDT/ParticleID.hh"


using namespace Rivet;


int WZandh::compare(const Projection & p) const {
  // they have to be same because there are no parameters.
  return true;
}

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
      HepPDT::ParticleID pInfo = (*pi)->pdg_id();
      int id = pInfo.abspid();
      if (id == 23 || id == 24 || id == 25 ){
	// This is a W Z or h.
	// Now find out whether it is the last one before the decay.
	// @todo trace the decay products of the bosons properly.
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
