// -*- C++ -*-
#ifndef RIVET_Particle_HH
#define RIVET_Particle_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/Vectors.hh"


namespace Rivet {
  
  /// Representation of particles from a HepMC::GenEvent.
  class Particle {
  public:    
    /// Default constructor.
    Particle() :
      _original(0), _id(0), _momentum(), _mass(0.0) 
    { }
    
    /// Constructor from a HepMC GenParticle.
    Particle(const GenParticle& gp) : 
      _original(&gp), _id(gp.pdg_id()), 
      _momentum(gp.momentum()), _mass(gp.momentum().m()) 
    { }

  public:
    /// Get a const reference to the original GenParticle.
    const GenParticle& getHepMCParticle() const { 
      assert(_original); 
      return *_original; 
    }
    
    /// Check if the particle corresponds to a GenParticle.
    bool hasHepMCParticle() const { return bool(_original); }
    
    /// The PDG ID code for this Particle.
    const long getPdgId() const { return _id; }

    /// The momentum of this Particle.
    FourMomentum& getMomentum() { return _momentum; }

    /// The momentum of this Particle.
    const FourMomentum& getMomentum() const { return _momentum; }

    /// Set the momentum of this Particle.
    Particle& setMomentum(const FourMomentum& momentum) { _momentum = momentum; return *this; }
    
    /// The mass of this Particle.
    const double getMass() const { return _mass; }

  
  private:
    /// A pointer to the original GenParticle from which this Particle is projected.
    const GenParticle* _original;
    
    /// The PDG ID code for this Particle.
    long _id;
    
    /// The momentum of this projection of the Particle.
    FourMomentum _momentum;
    
    /// @todo Check this reasoning!
    /// The mass of this Particle, stored for numerical hygiene.
    double _mass;    
  };
  
}

#endif
