// -*- C++ -*-
#ifndef RIVET_WZandh_HH
#define RIVET_WZandh_HH

#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {

  
  /// Project out the vector bosons and Higgs particles.
  /// @todo Remove this projection: it's unphysical and only used in one place (the ExampleTree)
  class WZandh : public Projection {
    
  public:
    
    /// @name Standard constructors and destructors. */
    //@{
    /// The default constructor. 
    WZandh() {
      setName("WZandh");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new WZandh(*this);
    }
    //@}
    
  protected:
    
    /// Apply the projection to the event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
    
  public:
    
    /// Access the Z's decayed to e+e-
    const ParticleVector& Zees() const { return _theZees; }

    /// Access the W's decayed to e neutrino
    const ParticleVector& Wens() const { return _theWens; }

    /// Access the W's decayed to mu neutrino
    const ParticleVector& Wmns() const { return _theWmns; }

    /// Access the Z's decayed to mu mu
    const ParticleVector& Zmms() const { return _theZmms; }

    /// Access the W's decayed to tau neutrino
    const ParticleVector& Wtns() const { return _theWtns; }

    /// Access the Z's decayed to tau tau
    const ParticleVector& Ztts() const { return _theZtts; }

    /// Access the Z's decayed to nu nu
    const ParticleVector& Znns() const { return _theZnns; }

    /// Access the W's decayed to qq
    const ParticleVector& Wqqs() const { return _theWqqs; }

    /// Access the Z's decayed to qq
    const ParticleVector& Zqqs() const { return _theZqqs; }

    /// Access Higgses
    const ParticleVector& hs() const { return _thehs; }
    
  private:
    
    /// @name The particles.
    //@{
    ParticleVector _theWens;
    ParticleVector _theZees;
    ParticleVector _theWmns;
    ParticleVector _theZmms;
    ParticleVector _theWtns;
    ParticleVector _theZtts;
    ParticleVector _theZnns;
    ParticleVector _theWqqs;
    ParticleVector _theZqqs;
    ParticleVector _thehs;
    //@}
        
  };
  
}


#endif
