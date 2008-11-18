// -*- C++ -*-
#ifndef RIVET_Particle_HH
#define RIVET_Particle_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/ParticleBase.hh"
#include "Rivet/Math/Vectors.hh"

namespace Rivet {

  
  /// Representation of particles from a HepMC::GenEvent.
  class Particle : public ParticleBase {
  public:    

    /// Default constructor.
    Particle() : ParticleBase(),
      _original(0), _id(0), _momentum(), _mass(0.0)
    { }

    
    /// Constructor from a HepMC GenParticle.
    Particle(const GenParticle& gp) : ParticleBase(),
    _original(&gp), _id(gp.pdg_id()),
    _momentum(gp.momentum()), _mass(gp.momentum().m())
    { }

  public:

    /// Get a const reference to the original GenParticle.
    const GenParticle& genParticle() const { 
      assert(_original); 
      return *_original; 
    }

    
    /// Check if the particle corresponds to a GenParticle.
    bool hasGenParticle() const { return bool(_original); }

    
    /// The PDG ID code for this Particle.
    const long pdgId() const { return _id; }


    /// The momentum of this Particle.
    FourMomentum& momentum() { return _momentum; }


    /// The momentum of this Particle.
    const FourMomentum& momentum() const { return _momentum; }


    /// Set the momentum of this Particle.
    Particle& setMomentum(const FourMomentum& momentum) { _momentum = momentum; return *this; }

    
    /// The mass of this Particle.
    const double mass() const { return _mass; }

    
    bool hasAncestor(long pdg_id) const {
      GenVertex* prodVtx = genParticle().production_vertex();
      if (prodVtx == 0) return false;
      GenVertex::particle_iterator ibegin = prodVtx->particles_begin(HepMC::ancestors);
      GenVertex::particle_iterator iend = prodVtx->particles_end(HepMC::ancestors);
      for (GenVertex::particle_iterator ancestor = ibegin; ancestor != iend; ++ancestor) {
        if ((*ancestor)->pdg_id() == pdg_id) {
          return true;
        }
      }      
      return false;
    }
    

  private:

    /// A pointer to the original GenParticle from which this Particle is projected.
    const GenParticle* _original;
    
    /// The PDG ID code for this Particle.
    long _id;
    
    /// The momentum of this projection of the Particle.
    FourMomentum _momentum;
    
    /// The mass of this Particle, stored for numerical hygiene.
    double _mass;    
  };

  
}

#endif
