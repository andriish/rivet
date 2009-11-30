// -*- C++ -*-
#ifndef RIVET_InitialQuarks_HH
#define RIVET_InitialQuarks_HH

#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// Project out quarks from the hard process in e+e- -> Z0 events
  /// @deprecated This is a very dangerous and specific projection! Use e.g. PID::hasBottom and friends instead
  class InitialQuarks : public Projection {
 
  public:
 
    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    InitialQuarks()
    {
      setName("InitialQuarks");
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new InitialQuarks(*this);
    }
    //@}
     
    /// Access the projected final-state particles.
    virtual const ParticleVector& particles() const { return _theParticles; }

    /// Is this final state empty?
    virtual const bool empty() const { return _theParticles.empty(); }

  protected:
 
    /// Apply the projection to the event.
    virtual void project(const Event& e);
 
    /// Compare projections.
    virtual int compare(const Projection& p) const;
 
  protected:

    /// The final-state particles.
    ParticleVector _theParticles;
 
  };

}


#endif
