// -*- C++ -*-
#include "Rivet/Particle.hh"
#include "Rivet/RivetCLHEP.hh"


namespace Rivet {

  Particle::Particle(const GenParticle& gp)
    : _original(&gp), _id(gp.pdg_id()), 
      _momentum(0), _mass(0.) {
    _momentum = new CLHEP::LorentzVector(gp.momentum().x(), gp.momentum().y(), 
                                         gp.momentum().z(), gp.momentum().t());
    _mass = gp.momentum().m();
  }
  
  
  Particle::Particle(const Particle& p)
    : _original(p._original), _id(p._id), 
      _momentum(0), _mass(p._mass) 
  { 
    if (p._momentum) {
      _momentum = new CLHEP::LorentzVector(*p._momentum);
    } else if (p._original) {
      const GenParticle& gp = *p._original;
      _momentum = new CLHEP::LorentzVector(gp.momentum().x(), gp.momentum().y(), 
                                           gp.momentum().z(), gp.momentum().t());
    }
  }
  
  
  Particle::~Particle() {
    delete _momentum;
  }


  Particle& Particle::operator=(const Particle& p) {
    if (&p != this) {
      _original = p._original;
      _id = p._id;
      _mass = p._mass;

      if (p._momentum) {
        _momentum = new CLHEP::LorentzVector(*p._momentum);
      } else if (p._original) {
        const GenParticle& gp = *p._original;
        _momentum = new CLHEP::LorentzVector(gp.momentum().x(), gp.momentum().y(), 
                                             gp.momentum().z(), gp.momentum().t());
      }
    }
    return *this;
  }


}
