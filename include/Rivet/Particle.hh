// -*- C++ -*-
#ifndef RIVET_Particle_HH
#define RIVET_Particle_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/ParticleBase.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {


  /// Representation of particles from a HepMC::GenEvent.
  class Particle : public ParticleBase {
  public:

    /// Default constructor.
    /// @deprecated A particle without info is useless. This only exists to keep STL containers happy.
    Particle() : ParticleBase(),
      _original(0), _id(0), _momentum()
    { }

    /// Constructor without GenParticle.
    Particle(PdgId pid, const FourMomentum& mom) : ParticleBase(),
      _original(0), _id(pid), _momentum(mom)
    { }

    /// Constructor from a HepMC GenParticle.
    Particle(const GenParticle& gp) : ParticleBase(),
    _original(&gp), _id(gp.pdg_id()),
    _momentum(gp.momentum())
    { }


  public:

    /// Get a const reference to the original GenParticle.
    const GenParticle& genParticle() const {
      assert(_original);
      return *_original;
    }


    /// Check if the particle corresponds to a GenParticle.
    bool hasGenParticle() const { 
      return bool(_original); 
    }


    /// The PDG ID code for this Particle.
    long pdgId() const { 
      return _id; 
    }


    /// The momentum of this Particle.
    const FourMomentum& momentum() const { 
      return _momentum; 
    }


    /// Set the momentum of this Particle.
    Particle& setMomentum(const FourMomentum& momentum) { 
      _momentum = momentum; 
      return *this; 
    }


    /// The mass of this Particle.
    double mass() const { 
      return momentum().mass(); 
    }

    // /// The charge of this Particle.
    // double charge() const { 
    //   return PID::charge(*this); 
    // }

    // /// Three times the charge of this Particle (i.e. integer multiple of smallest quark charge).
    // int threeCharge() const { 
    //   return PID::threeCharge(*this); 
    // }


    /// Check whether a given PID is found in the GenParticle's ancestor list
    bool hasAncestor(PdgId pdg_id) const;


  private:

    /// A pointer to the original GenParticle from which this Particle is projected.
    const GenParticle* _original;

    /// The PDG ID code for this Particle.
    long _id;

    /// The momentum of this projection of the Particle.
    FourMomentum _momentum;
  };



  inline bool cmpParticleByPt(const Particle& a, const Particle& b) {
    return a.momentum().pT() > b.momentum().pT();
  }
  inline bool cmpParticleByAscPt(const Particle& a, const Particle& b) {
    return a.momentum().pT() < b.momentum().pT();
  }
  inline bool cmpParticleByEt(const Particle& a, const Particle& b) {
    return a.momentum().Et() > b.momentum().Et();
  }
  inline bool cmpParticleByAscEt(const Particle& a, const Particle& b) {
    return a.momentum().Et() < b.momentum().Et();
  }
  inline bool cmpParticleByE(const Particle& a, const Particle& b) {
    return a.momentum().E() > b.momentum().E();
  }
  inline bool cmpParticleByAscE(const Particle& a, const Particle& b) {
    return a.momentum().E() < b.momentum().E();
  }
  inline bool cmpParticleByDescPseudorapidity(const Particle& a, const Particle& b) {
    return a.momentum().pseudorapidity() > b.momentum().pseudorapidity();
  }
  inline bool cmpParticleByAscPseudorapidity(const Particle& a, const Particle& b) {
    return a.momentum().pseudorapidity() < b.momentum().pseudorapidity();
  }
  inline bool cmpParticleByDescAbsPseudorapidity(const Particle& a, const Particle& b) {
    return fabs(a.momentum().pseudorapidity()) > fabs(b.momentum().pseudorapidity());
  }
  inline bool cmpParticleByAscAbsPseudorapidity(const Particle& a, const Particle& b) {
    return fabs(a.momentum().pseudorapidity()) < fabs(b.momentum().pseudorapidity());
  }
  inline bool cmpParticleByDescRapidity(const Particle& a, const Particle& b) {
    return a.momentum().rapidity() > b.momentum().rapidity();
  }
  inline bool cmpParticleByAscRapidity(const Particle& a, const Particle& b) {
    return a.momentum().rapidity() < b.momentum().rapidity();
  }
  inline bool cmpParticleByDescAbsRapidity(const Particle& a, const Particle& b) {
    return fabs(a.momentum().rapidity()) > fabs(b.momentum().rapidity());
  }
  inline bool cmpParticleByAscAbsRapidity(const Particle& a, const Particle& b) {
    return fabs(a.momentum().rapidity()) < fabs(b.momentum().rapidity());
  }




}

#endif
