// -*- C++ -*-
#ifndef RIVET_Particle_HH
#define RIVET_Particle_HH

#include "Rivet/Rivet.hh"
#include "HepMC/GenParticle.h"
#include "Rivet/Particle.fhh"
#include "Rivet/RivetCLHEP.fhh"


namespace Rivet {
  /// Typedef from HepMC.
  typedef HepMC::GenParticle GenParticle;

  
  /// This class is used to represent particles projected from a HepMC::GenEvent.
  class Particle {
  public:

    /// Default constructor.
    inline Particle() :
      _original(0), _id(0), _momentum(0), _mass(0.0) 
    { };
    
    /// Constructor from a HepMC GenParticle.
    Particle(const GenParticle& gp);
    
    /// Copy-constructor.
    Particle(const Particle& p);

    /// Destructor.
    ~Particle();

    // Copy-assignment.
    Particle& operator=(const Particle& p);


  public:
    /// Get a const reference to the original GenParticle.
    inline const GenParticle& getHepMCParticle() const { return *_original; }
    
    /// Check if the particle corresponds to a GenParticle.
    inline bool hasHepMCParticle() const { return _original; }
    
    /// The PDG ID code for this Particle.
    inline const long getPdgId() const { return _id; }
    
    /// The momentum of this projection of the Particle.
    inline CLHEP::LorentzVector& getMomentum() { return *_momentum; }

    /// The momentum of this projection of the Particle.
    inline const CLHEP::LorentzVector& getMomentum() const { return *_momentum; }
    
    /// The mass of this Particle.
    inline const double getMass() const { return _mass; }

  
  private:
  
    /// A pointer to the original GenParticle from which this Particle is projected.
    const GenParticle* _original;
    
    /// The PDG ID code for this Particle.
    long _id;
    
    /// The momentum of this projection of the Particle.
    CLHEP::LorentzVector* _momentum;
    
    /// The mass of this Particle.
    double _mass;
    
  };
  
}

#endif
